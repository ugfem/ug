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
#include "general.h"
#include "ugm.h"      /* for GetFreeOBJT() */
#include "memmgr.h"

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

/* macro for easy definition of type mapping UG<->DDD */
#define MAP_TYPES(ugt,dddt)   { int _ugt=(ugt); \
                                dddctrl.ugtypes[(dddt)] = _ugt;     \
                                dddctrl.types[_ugt] = (dddt);       \
}


/* #define DDD_PrioMergeDefault(x,y)*/	/* TODO: delete this define */


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
#ifdef __THREEDIM__
DDD_TYPE TypeEdge;
#endif

DDD_TYPE TypeUnknown;

#ifdef __TWODIM__
DDD_TYPE TypeTrElem, TypeTrBElem;
DDD_TYPE TypeQuElem, TypeQuBElem;
#endif

#ifdef __THREEDIM__
DDD_TYPE TypeTeElem, TypeTeBElem;
DDD_TYPE TypePyElem, TypePyBElem;
DDD_TYPE TypePrElem, TypePrBElem;
DDD_TYPE TypeHeElem, TypeHeBElem;
#endif


/* DDD data objects */
DDD_TYPE TypeMatrix;
DDD_TYPE TypeBndP;
#ifdef __TWODIM__
DDD_TYPE TypeEdge;
#endif
DDD_TYPE TypeBndS;

/* DDD interfaces needed for distributed computation */
DDD_IF ElementIF, ElementSymmIF, ElementVIF, ElementSymmVIF,
       ElementVHIF, ElementSymmVHIF;
DDD_IF BorderNodeIF, BorderNodeSymmIF, OuterNodeIF, NodeVIF,
       NodeIF, NodeAllIF;
DDD_IF BorderVectorIF, BorderVectorSymmIF,
       OuterVectorIF, OuterVectorSymmIF,
       VectorVIF, VectorVAllIF, VectorAllIF;
DDD_IF VertexIF;
#ifdef __THREEDIM__
DDD_IF EdgeIF, BorderEdgeSymmIF, EdgeHIF, EdgeAllIF;
#endif



/* DDD global controls */
DDD_CTRL dddctrl;


/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

enum ElemTypeFlag { Inside, Boundary };


/****************************************************************************/
/*
   void ddd_InitGenericElement -

   SYNOPSIS:
   static void ddd_InitGenericElement (INT tag, DDD_TYPE dddType, int etype);

   PARAMETERS:
   .  tag
   .  dddType
   .  etype

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void ddd_InitGenericElement (INT tag, DDD_TYPE dddType, int etype)
{
  struct generic_element  *ge=0;
  GENERAL_ELEMENT                 *desc = element_descriptors[tag];
  unsigned INT gbits = 0;

  size_t ps  = sizeof(void *);
  void   **r = ge->refs;

  /* compute global fields of control word entry */
  gbits = ~(((1<<NSONS_LEN)-1)<<NSONS_SHIFT);
  PRINTDEBUG(dddif,1,("ddd_InitGenericElement(): gbits=%08x size=%d\n",
                      gbits,sizeof(ge->control)));

  /* initialize base part (valid for all elements) */
  DDD_TypeDefine(dddType, ge,
                 EL_DDDHDR, &(ge->ddd),
                 /* TODO: delete this					*/
                 /*		EL_GDATA,  ELDEF(ge->control),	*/
                 EL_GBITS,  ELDEF(ge->control), &gbits,

                 /* TODO: id muss umgerechnet werden! (?) */
                 EL_GDATA,  ELDEF(ge->id),
                 EL_GDATA,  ELDEF(ge->flag),
                 EL_GDATA,  ELDEF(ge->property),
                 EL_GDATA,  ELDEF(ge->ptmp),
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
                 EL_OBJPTR, r+father_offset[tag],  ps,                       dddType,
                 /* TODO: delete
                    #ifdef __TWODIM__
                                 EL_LDATA, r+sons_offset[tag],    ps*desc->max_sons_of_elem,
                    #endif
                    #ifdef __THREEDIM__
                                 EL_LDATA, r+sons_offset[tag],    ps*1,
                    #endif
                  */
                 EL_LDATA, r+sons_offset[tag],    ps*2,
                 EL_OBJPTR, r+nb_offset[tag],      ps*desc->sides_of_elem,   dddType,
                 EL_CONTINUE);


  /* optional components */

  if (dddctrl.elemData)
    DDD_TypeDefine(dddType, ge,
                   EL_OBJPTR, r+evector_offset[tag], ps*1,     TypeVector,
                   EL_CONTINUE);

        #ifdef __THREEDIM__
  if (dddctrl.sideData)
    DDD_TypeDefine(dddType, ge,
                   EL_OBJPTR, r+svector_offset[tag], ps*desc->sides_of_elem, TypeVector,
                   EL_CONTINUE);
        #endif

  if (etype==Inside)
  {
    DDD_TypeDefine(dddType, ge, EL_END, desc->inner_size);

    /* init type mapping arrays */
    MAP_TYPES(MAPPED_INNER_OBJT_TAG(tag), dddType);
    dddctrl.dddObj[MAPPED_INNER_OBJT_TAG(tag)] = TRUE;
  }
  else
  {
    DDD_TypeDefine(dddType, ge,
                   EL_LDATA, r+side_offset[tag],  ps*desc->sides_of_elem,
                   EL_END, desc->bnd_size);

    /* init type mapping arrays */
    MAP_TYPES(MAPPED_BND_OBJT_TAG(tag), dddType);
    dddctrl.dddObj[MAPPED_BND_OBJT_TAG(tag)] = TRUE;
  }

  /* set mergemode to maximum */
  DDD_PrioMergeDefault(dddType, PRIOMERGE_MAXIMUM);
  /* TODO: set prios
          DDD_PrioMergeDefine(dddType, PrioHGhost, PrioVGhost, PrioVHGhost);
          DDD_PrioMergeDefine(dddType, PrioHGhost, PrioVHGhost, PrioVHGhost);
          DDD_PrioMergeDefine(dddType, PrioVGhost, PrioVHGhost, PrioVHGhost);
          DDD_PrioMergeDefine(dddType, PrioHGhost, PrioMaster, PrioMaster);
          DDD_PrioMergeDefine(dddType, PrioVGhost, PrioMaster, PrioMaster);
          DDD_PrioMergeDefine(dddType, PrioVHGhost, PrioMaster, PrioMaster);
   */

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


/****************************************************************************/
/*
   void ddd_DeclareTypes - declare ug data structures as DDD_TYPES

   SYNOPSIS:
   static void ddd_DeclareTypes (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function declares ug data structures as DDD_TYPES

   RETURN VALUE:
   void
 */
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
  MAP_TYPES(VEOBJ, TypeVector);
  dddctrl.dddObj[VEOBJ] = TRUE;

  TypeIVertex     = DDD_TypeDeclare("IVertex");
  MAP_TYPES(IVOBJ, TypeIVertex);
  dddctrl.dddObj[IVOBJ] = TRUE;

  TypeBVertex     = DDD_TypeDeclare("BVertex");
  MAP_TYPES(BVOBJ, TypeBVertex);
  dddctrl.dddObj[BVOBJ] = TRUE;

  TypeNode        = DDD_TypeDeclare("Node");
  MAP_TYPES(NDOBJ, TypeNode);
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
  TypePrElem      = DDD_TypeDeclare("PrElem");
  TypePrBElem     = DDD_TypeDeclare("PrBElem");
  TypeHeElem      = DDD_TypeDeclare("HeElem");
  TypeHeBElem     = DDD_TypeDeclare("HeBElem");
        #endif /* THREEDIM */

  /* edge type not unique:                    */
  /* edge is DDD object for 3D                */
  /* edge is DDD data object for 2D           */
  TypeEdge        = DDD_TypeDeclare("Edge");
  MAP_TYPES(EDOBJ, TypeEdge);
        #ifdef __THREEDIM__
  dddctrl.dddObj[EDOBJ] = TRUE;
        #endif

  /* 2. DDD data objects (without DDD_HEADER) */

  TypeMatrix  = DDD_TypeDeclare("Matrix");
  MAP_TYPES(MAOBJ, TypeMatrix);

  TypeBndP    = DDD_TypeDeclare("BndP");
  MAP_TYPES(GetFreeOBJT(), TypeBndP);

  TypeBndS = DDD_TypeDeclare("BndS");
  MAP_TYPES(GetFreeOBJT(), TypeBndS);
}


/****************************************************************************/
/*
   ddd_DefineTypes - define previously declared DDD_TYPES

   SYNOPSIS:
   static void ddd_DefineTypes (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function defines previously declared DDD_TYPES
   Note: this function depends on previous definition of all necessary ug-generic-elements.

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void ddd_DefineTypes (void)
{
  INT size;
  VECTOR v;
  NODE n;
  struct ivertex iv;
  struct bvertex bv;

  MATRIX m;
  EDGE e;
  unsigned INT gbits = 0;

  /* 1. DDD objects (with DDD_HEADER) */

  DDD_TypeDefine(TypeVector, &v,
                 EL_DDDHDR, &v.ddd,
                 EL_GDATA,  ELDEF(v.control),

                 /* object must be LDATA, because reftype may be a non-DDD-object */
                 /* (e.g., edge). therefore, 'object' must be updated by MKCONS-  */
                 /* handler of associated object. 960404 KB */
                 /* TODO: decide whether LDATA or OBJPTR for different VectorTypes*/
                 /* EL_OBJPTR, ELDEF(v.object), TypeNode, */
                 EL_LDATA,  ELDEF(v.object),
                 EL_LDATA,  ELDEF(v.pred),
                 EL_LDATA,  ELDEF(v.succ),
                 EL_GDATA,  ELDEF(v.index),
                 EL_GDATA,  ELDEF(v.skip),
                 EL_LDATA,  ELDEF(v.start),

                 /* TODO is this LDATA or GDATA? */
                #ifdef __BLOCK_VECTOR_DESC__
                 EL_LDATA,  ELDEF(v.block_descr),
                #endif

                #ifdef __INTERPOLATION_MATRIX__
                 EL_LDATA,  ELDEF(v.istart),
                #endif

                 /* TODO: value wird noch ausgelassen. feld variabler laenge? */
                 /* bei entscheidung 'value': kein weiteres feld
                         bei ent. 'userdata *': EL_GDATA-feld        */
                 EL_GDATA,  ELDEF(v.value),
                 EL_END,    &v+1
                 );

  /* set mergemode to maximum */
  DDD_PrioMergeDefault(TypeVector, PRIOMERGE_MAXIMUM);

  /* compute global fields it control word entry */
  gbits = ~((((1<<ONEDGE_LEN)-1)<<ONEDGE_SHIFT) |
            (((1<<NOOFNODE_LEN)-1)<<NOOFNODE_SHIFT));
  PRINTDEBUG(dddif,1,("ddd_DefineTypes(): TypeI/BVertex gbits=%08x size=%d\n",
                      gbits,sizeof(iv.control)));

  DDD_TypeDefine(TypeIVertex, &iv,
                 EL_DDDHDR, &iv.ddd,
                 /* TODO: delete
                                 EL_GDATA,  ELDEF(iv.control),
                  */
                 EL_GBITS,  ELDEF(iv.control), &gbits,

                 /* TODO: muss umgerechnet werden! */
                 EL_GDATA,  ELDEF(iv.id),
                 EL_GDATA,  ELDEF(iv.x),
                 EL_GDATA,  ELDEF(iv.xi),
                 EL_LDATA,  ELDEF(iv.pred),
                 EL_LDATA,  ELDEF(iv.succ),
                 EL_LDATA,  ELDEF(iv.data),

                 /* TODO muss father LDATA oder OBJPTR sein?     */
                 /* LDATA, father ist nur lokal gueltig und      */
                 /* ist abhaengig von vertikaler Lastverteilung  */
                #ifdef __TWODIM__
                 /* TODO: ref-typ muss eigentlich {TypeTrElem,TypeTrBElem} sein! */
                 EL_OBJPTR, ELDEF(iv.father), TypeTrElem,
                #endif
                #ifdef __THREEDIM__
                 EL_LDATA,  ELDEF(iv.father),
                #endif

                #ifdef TOPNODE
                 /* TODO topnode wirklich OBJPTR? */
                 EL_LDATA,  ELDEF(iv.topnode),
                #endif
                 EL_END,    &iv+1
                 );

  /* set mergemode to maximum */
  DDD_PrioMergeDefault(TypeIVertex, PRIOMERGE_MAXIMUM);


  DDD_TypeDefine(TypeBVertex, &bv,
                 EL_DDDHDR, &bv.ddd,
                 /* TODO: delete
                                 EL_GDATA,  ELDEF(bv.control),
                  */
                 EL_GBITS,  ELDEF(bv.control), &gbits,

                 /* TODO: muss umgerechnet werden! Nooeee! */
                 EL_GDATA,  ELDEF(bv.id),
                 EL_GDATA,  ELDEF(bv.x),
                 EL_GDATA,  ELDEF(bv.xi),
                 EL_LDATA,  ELDEF(bv.pred),
                 EL_LDATA,  ELDEF(bv.succ),
                 EL_LDATA,  ELDEF(bv.data),

                 /* TODO muss father LDATA oder OBJPTR sein?     */
                 /* LDATA, father ist nur lokal gueltig und      */
                 /* ist abhaengig von vertikaler Lastverteilung  */
                #ifdef __TWODIM__
                 /* TODO: ref-typ muss eigentlich {TypeTrElem,TypeTrBElem} sein! */
                 EL_OBJPTR, ELDEF(bv.father), TypeTrElem,
                #endif
                #ifdef __THREEDIM__
                 EL_LDATA,  ELDEF(bv.father),
                #endif

                #ifdef TOPNODE
                 /* TODO topnode wirklich OBJPTR?, Nooeee! */
                 EL_LDATA,  ELDEF(bv.topnode),
                #endif
                 EL_LDATA,  ELDEF(bv.bndp),     /* different from IVertex */
                 EL_END,    &bv+1
                 );

  /* set mergemode to maximum */
  DDD_PrioMergeDefault(TypeBVertex, PRIOMERGE_MAXIMUM);


  DDD_TypeDefine(TypeNode, &n,
                 EL_DDDHDR, &n.ddd,
                 EL_GDATA,  ELDEF(n.control),

                 /* TODO: muss umgerechnet werden! */
                 EL_GDATA,  ELDEF(n.id),
                 EL_LDATA,  ELDEF(n.pred),
                 EL_LDATA,  ELDEF(n.succ),

                 /* TODO was ist start? */
                 EL_LDATA,  ELDEF(n.start),

                 /* father may be one of node or edge */
                 EL_OBJPTR, ELDEF(n.father),   DDD_TYPE_BY_HANDLER, NFatherObjType,

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

  /* set mergemode to maximum */
  DDD_PrioMergeDefault(TypeNode, PRIOMERGE_MAXIMUM);


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
  ddd_InitGenericElement(PRISM,           TypePrElem,  Inside);
  ddd_InitGenericElement(PRISM,           TypePrBElem, Boundary);
  ddd_InitGenericElement(HEXAHEDRON,  TypeHeElem,  Inside);
  ddd_InitGenericElement(HEXAHEDRON,  TypeHeBElem, Boundary);
        #endif /* THREEDIM */

  /* 2. DDD data objects (without DDD_HEADER) */

  /* NOTE: the size of matrix objects computed by the DDD Typemanager
     will not be the real size ued by DDD. this size has to be computed
     by MSIZE(mat). this is relevant only in gather/scatter of matrices
     in handler.c. */
  DDD_TypeDefine(TypeMatrix, &m,
                 EL_GDATA,  ELDEF(m.control),
                 EL_LDATA,  ELDEF(m.next),
                 EL_OBJPTR, ELDEF(m.vect),   TypeVector,
                 /* TODO: not needed
                    EL_LDATA,  ELDEF(m.value), */
                 EL_END,    &m+1
                 );

  /* compute global fields it control word entry */
  gbits = ~(((1<<NO_OF_ELEM_LEN)-1)<<NO_OF_ELEM_SHIFT);
  PRINTDEBUG(dddif,1,("ddd_DefineTypes(): TypeEdge gbits=%08x size=%d\n",
                      gbits,sizeof(e.links[0].control)));

  DDD_TypeDefine(TypeEdge, &e,
                 /* link 0 data */
                 /*TODO: now unique
                    #ifdef __TWODIM__
                                 EL_GDATA,  ELDEF(e.links[0].control),
                    #endif
                    #ifdef __THREEDIM__
                                 EL_LDATA,  ELDEF(e.links[0].control),
                    #endif
                  */
                 EL_GBITS,  ELDEF(e.links[0].control), &gbits,
                 EL_LDATA,  ELDEF(e.links[0].next),
                 EL_OBJPTR, ELDEF(e.links[0].nbnode), TypeNode,

                 /* link 1 data */
                 EL_GDATA,  ELDEF(e.links[1].control),
                 EL_LDATA,  ELDEF(e.links[1].next),
                 EL_OBJPTR, ELDEF(e.links[1].nbnode), TypeNode,

                #ifdef __THREEDIM__
                 EL_DDDHDR, &e.ddd,
                #endif

                 EL_OBJPTR, ELDEF(e.midnode),  TypeNode,
                 EL_CONTINUE);

  if (dddctrl.edgeData)
    DDD_TypeDefine(TypeEdge, &e,
                   EL_OBJPTR, ELDEF(e.vector), TypeVector,
                   EL_CONTINUE);

  size = sizeof(EDGE) - ((dddctrl.edgeData) ? 0 : sizeof(VECTOR*));
  DDD_TypeDefine(TypeEdge, &e, EL_END, ((char *)&e)+size);

        #ifdef __THREEDIM__
  /* set mergemode to maximum */
  DDD_PrioMergeDefault(TypeEdge, PRIOMERGE_MAXIMUM);
        #endif

}


/****************************************************************************/
/*
   ddd_IfInit - define the communication interfaces needed in ug for management by DDD

   SYNOPSIS:
   static void ddd_IfInit (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function defines the communication interfaces needed in ug for management by DDD

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void ddd_IfInit (void)
{
  DDD_TYPE O[8];
  int nO;
  DDD_PRIO A[8];
  DDD_PRIO B[8];


  /* define element interfaces */
#ifdef __TWODIM__
  O[0] = TypeTrElem; O[1] = TypeTrBElem;
  O[2] = TypeQuElem; O[3] = TypeQuBElem;
  nO = 4;
#endif

#ifdef __THREEDIM__
  O[0] = TypeTeElem; O[1] = TypeTeBElem;
  O[2] = TypePyElem; O[3] = TypePyBElem;
  O[4] = TypePrElem; O[5] = TypePrBElem;
  O[6] = TypeHeElem; O[7] = TypeHeBElem;
  nO = 8;
#endif

  A[0] = PrioMaster;
  B[0] = PrioHGhost; B[1] = PrioVHGhost;
  ElementIF = DDD_IFDefine(nO,O,1,A,2,B);
  DDD_IFSetName(ElementIF, "ElementIF: Master->HGhost");

  A[0] = PrioMaster; A[1] = PrioHGhost; A[2] = PrioVHGhost;
  B[0] = PrioMaster; B[1] = PrioHGhost; B[2] = PrioVHGhost;
  ElementSymmIF = DDD_IFDefine(nO,O,3,A,3,B);
  DDD_IFSetName(ElementSymmIF, "ElementSymmIF: Master/HGhost");

  A[0] = PrioMaster;
  B[0] = PrioVGhost; B[1] = PrioVHGhost;
  ElementVIF = DDD_IFDefine(nO,O,1,A,2,B);
  DDD_IFSetName(ElementVIF, "ElementVIF: Master->VGhost");

  A[0] = PrioMaster; A[1] = PrioVGhost; A[2] = PrioVHGhost;
  B[0] = PrioMaster; B[1] = PrioVGhost; B[2] = PrioVHGhost;
  ElementSymmVIF = DDD_IFDefine(nO,O,3,A,3,B);
  DDD_IFSetName(ElementSymmVIF, "ElementSymmVIF: Master/VGhost");

  A[0] = PrioMaster;
  B[0] = PrioVGhost; B[1] = PrioHGhost; B[2] = PrioVHGhost;
  ElementVHIF = DDD_IFDefine(nO,O,1,A,3,B);
  DDD_IFSetName(ElementVHIF, "ElementVHIF: Master->VGhost/HGhost");

  A[0] = PrioMaster; A[1] = PrioVGhost; A[2] = PrioHGhost; A[3] = PrioVHGhost;
  B[0] = PrioMaster; B[1] = PrioVGhost; B[2] = PrioHGhost; B[3] = PrioVHGhost;
  ElementSymmVHIF = DDD_IFDefine(nO,O,4,A,4,B);
  DDD_IFSetName(ElementSymmVHIF, "ElementSymmVHIF: Master/VGhost/HGhost");


  /* define node interfaces */
  O[0] = TypeNode;

  A[0] = PrioBorder;
  B[0] = PrioMaster;
  BorderNodeIF = DDD_IFDefine(1,O,1,A,1,B);
  DDD_IFSetName(BorderNodeIF, "BorderNodeIF: Border->Master");

  A[0] = PrioMaster; A[1] = PrioBorder;
  B[0] = PrioMaster; B[1] = PrioBorder;
  BorderNodeSymmIF = DDD_IFDefine(1,O,2,A,2,B);
  DDD_IFSetName(BorderNodeSymmIF, "BorderNodeSymmIF: Border/Master");

  A[0] = PrioMaster;
  B[0] = PrioHGhost; B[1] = PrioVHGhost;
  OuterNodeIF = DDD_IFDefine(1,O,1,A,2,B);
  DDD_IFSetName(OuterNodeIF, "OuterNodeIF: Master->HGhost");

  A[0] = PrioMaster;
  B[0] = PrioVGhost; B[1] = PrioVHGhost;
  NodeVIF = DDD_IFDefine(1,O,1,A,2,B);
  DDD_IFSetName(NodeVIF, "NodeVIF: Master->VGhost");

  A[0] = PrioMaster;
  B[0] = PrioVGhost; B[1] = PrioHGhost; B[2] = PrioVHGhost;
  NodeIF = DDD_IFDefine(1,O,1,A,3,B);
  DDD_IFSetName(NodeIF, "NodeIF: Master->HGhost");

  A[0] = PrioMaster; A[1] = PrioBorder; A[2] = PrioVGhost; A[3] = PrioHGhost; A[4] = PrioVHGhost;
  B[0] = PrioMaster; B[1] = PrioBorder; B[2] = PrioVGhost; B[3] = PrioHGhost; B[4] = PrioVHGhost;
  NodeAllIF = DDD_IFDefine(1,O,5,A,5,B);
  DDD_IFSetName(NodeAllIF, "NodeAllIF: All/All");


  /* define vector interfaces */
  O[0] = TypeVector;

  A[0] = PrioBorder;
  B[0] = PrioMaster;
  BorderVectorIF = DDD_IFDefine(1,O,1,A,1,B);
  DDD_IFSetName(BorderVectorIF, "BorderVectorIF: Border->Master");

  A[0] = PrioMaster; A[1] = PrioBorder;
  B[0] = PrioMaster; B[1] = PrioBorder;
  BorderVectorSymmIF = DDD_IFDefine(1,O,2,A,2,B);
  DDD_IFSetName(BorderVectorSymmIF, "BorderVectorSymmIF: Master/Border");

  A[0] = PrioMaster;
  B[0] = PrioHGhost; B[1] = PrioVHGhost;
  OuterVectorIF = DDD_IFDefine(1,O,1,A,2,B);
  DDD_IFSetName(OuterVectorIF, "OuterVectorIF: Master->HGhost");

  A[0] = PrioMaster; A[1] = PrioBorder; A[2] = PrioHGhost; A[3] = PrioVHGhost;
  B[0] = PrioMaster; B[1] = PrioBorder; B[2] = PrioHGhost; B[3] = PrioVHGhost;
  OuterVectorSymmIF = DDD_IFDefine(1,O,4,A,4,B);
  DDD_IFSetName(OuterVectorSymmIF, "OuterVectorSymmIF: Master/Border/HGhost");

  A[0] = PrioMaster;
  B[0] = PrioVGhost; B[1] = PrioVHGhost;
  VectorVIF = DDD_IFDefine(1,O,1,A,2,B);
  DDD_IFSetName(VectorVIF, "VectorVIF: Master->VGhost");

  A[0] = PrioMaster; A[1] = PrioBorder; A[2] = PrioVGhost; A[3] = PrioVHGhost;
  B[0] = PrioMaster; B[1] = PrioBorder;
  VectorVAllIF = DDD_IFDefine(1,O,4,A,2,B);
  DDD_IFSetName(VectorVAllIF, "VectorVAllIF: Master/Border/VGhost->Master/Border");

  A[0] = PrioMaster;
  B[0] = PrioBorder; B[1] = PrioVGhost;
  B[2] = PrioVHGhost; B[3] = PrioHGhost;
  VectorAllIF = DDD_IFDefine(1,O,1,A,4,B);
  DDD_IFSetName(VectorVAllIF, "VectorAllIF: Master->Border/HGhost/VGhost/VHGhost");

  /* define vertex interfaces */
  O[0] = TypeIVertex; O[1] = TypeBVertex;

  A[0] = PrioMaster;
  B[0] = PrioMaster;
  VertexIF = DDD_IFDefine(2,O,1,A,1,B);
  DDD_IFSetName(VertexIF, "VertexIF: Master<->Master");


  /* define edge interfaces */
        #ifdef __THREEDIM__
  O[0] = TypeEdge;

  A[0] = PrioMaster;
  B[0] = PrioMaster;
  EdgeIF = DDD_IFDefine(1,O,1,A,1,B);
  DDD_IFSetName(EdgeIF, "EdgeIF: Master<->Master");

  A[0] = PrioMaster; A[1] = PrioBorder;
  B[0] = PrioMaster; B[1] = PrioBorder;
  BorderEdgeSymmIF = DDD_IFDefine(1,O,2,A,2,B);
  DDD_IFSetName(BorderEdgeSymmIF, "BorderEdgeSymmIF: Master/Border");

  A[0] = PrioMaster; A[1] = PrioBorder;
  B[0] = PrioMaster; B[1] = PrioBorder; B[2] = PrioHGhost; B[3] = PrioVHGhost;
  EdgeHIF = DDD_IFDefine(1,O,2,A,4,B);
  DDD_IFSetName(EdgeHIF, "EdgeHIF: Master/Border");

  A[0] = PrioMaster; A[1] = PrioBorder; A[2] = PrioVGhost; A[3] = PrioHGhost; A[4] = PrioVHGhost;
  B[0] = PrioMaster; B[1] = PrioBorder; B[2] = PrioVGhost; B[3] = PrioHGhost; B[4] = PrioVHGhost;
  EdgeAllIF = DDD_IFDefine(1,O,5,A,5,B);
  DDD_IFSetName(EdgeAllIF, "EdgeAllIF: All/All");

        #endif
}


/****************************************************************************/
/*
   InitDDDTypes - define DDD_TYPEs

   SYNOPSIS:
   static void InitDDDTypes (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function must be called once before creation of DDD-objects. It depends on correct and complete initialization of all ug-generic-elements, therefore it must be called after completion of InitElementTypes(). As InitElementTypes() will be called whenever new Multigrids are created/opened, an execution guard prevents this function from multiple execution.

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void InitDDDTypes (void)
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
        #ifdef __THREEDIM__
  DDD_TypeDisplay(TypeEdge);
        #endif

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
  DDD_TypeDisplay(TypePrElem);
  DDD_TypeDisplay(TypePrBElem);
  DDD_TypeDisplay(TypeHeElem);
  DDD_TypeDisplay(TypeHeBElem);
        #endif

  /* display dependent types */
  DDD_TypeDisplay(TypeMatrix);
        #ifdef __TWODIM__
  DDD_TypeDisplay(TypeEdge);
        #endif
  ENDDEBUG

  ddd_HandlerInit(HSET_XFER);
}



/****************************************************************************/
/*
   InitCurrMG - initialize the current multigrid which is handled by DDD

   SYNOPSIS:
   void InitCurrMG (MULTIGRID *MG);

   PARAMETERS:
   .  MG

   DESCRIPTION:
   This function initializes the current multigrid which is handled by DDD

   RETURN VALUE:
   void
 */
/****************************************************************************/

void InitCurrMG (MULTIGRID *MG)
{
  dddctrl.currMG = MG;

  dddctrl.nodeData = VEC_DEF_IN_OBJ_OF_MG(dddctrl.currMG,NODEVEC);
  dddctrl.edgeData = VEC_DEF_IN_OBJ_OF_MG(dddctrl.currMG,EDGEVEC);
  dddctrl.elemData = VEC_DEF_IN_OBJ_OF_MG(dddctrl.currMG,ELEMVEC);
  dddctrl.sideData = VEC_DEF_IN_OBJ_OF_MG(dddctrl.currMG,SIDEVEC);

  if (dddctrl.currFormat == NULL) {
    InitDDDTypes();
    dddctrl.currFormat = MGFORMAT(MG);
  }
  else if (dddctrl.currFormat !=  MGFORMAT(MG))
    assert(0);
}

/****************************************************************************/
/*
   CheckInitParallel - check for correct initialization of dddif subsystem

   SYNOPSIS:
   static int CheckInitParallel (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function performs checks for correct initialization of dddif subsystem

   RETURN VALUE:
   int
     error value
 */
/****************************************************************************/


static int CheckInitParallel (void)
{
  int i;

  /* check for valid UGTYPE for given DDD_TYPE */
  if (OBJT_MAX == MAXOBJECTS)
  {
    printf("ERROR in InitParallel: OBJT_MAX!=MAXOBJECTS\n");
    return(__LINE__);
  }

  for(i=1; i<MAXDDDTYPES && UGTYPE(i)>=0; i++)
  {
    /* check for valid UGTYPE for given DDD_TYPE */
    if (UGTYPE(i) > OBJT_MAX)
    {
      printf("ERROR in InitParallel: OBJT=%d > OBJT_MAX=%d\n",
             UGTYPE(i), OBJT_MAX);
      return(__LINE__);
    }

    /* check for correct mapping and re-mapping */
    if (DDDTYPE(UGTYPE(i))!=i)
    {
      printf("ERROR in InitParallel: invalid type mapping for OBJT=%d\n",
             UGTYPE(i));
      return(__LINE__);
    }
  }

  /* no errors */
  return(0);
}



/****************************************************************************/
/*
   InitParallel - initialize the ddd library

   SYNOPSIS:
   int InitParallel (int *argc, char ***argv);

   PARAMETERS:
   .  argc - pointer to number of arguments
   .  argv - pointer to list of argument pointers

   DESCRIPTION:
   This function initializes the ddd library by defining the ug internal issues: format of handled structs, description of handlers, definition of interfaces

   RETURN VALUE:
   int
 */
/****************************************************************************/

int InitParallel (void)
{
  INT err;
  int i;

  memmgr_Init();

  /* init DDD and set options */
  DDD_Init(NULL,NULL);

  /* we are using varsized DDD objects, turn warnings off */
  DDD_SetOption(OPT_WARNING_VARSIZE_OBJ, OPT_OFF);
  DDD_SetOption(OPT_WARNING_SMALLSIZE, OPT_OFF);

  /* show messages during transfer, for debugging */
  DDD_SetOption(OPT_DEBUG_XFERMESGS, OPT_OFF);

  /* TODO: remove this, reference collision with Edge orientation
     in 3D */
  DDD_SetOption(OPT_WARNING_REF_COLLISION, OPT_OFF);

  /* treat identify tokens for one object as set */
  DDD_SetOption(OPT_IDENTIFY_MODE, IDMODE_SETS);

  /* dont delete objects when another copy comes in during Xfer */
  DDD_SetOption(OPT_XFER_PRUNE_DELETE, OPT_ON);

  /* initialize context */
  /* TODO: malloc() should be replaced by HEAPs or ddd_memmgr */
  dddctrl._context = (INT *)malloc(sizeof(INT)*procs);

  /* initial context is all processors */
  for(i=0; i<procs; i++)
    dddctrl._context[i] = 1;


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
  dddctrl.currFormat = NULL;

  /* declare DDD_TYPES, definition must be done later */
  ddd_DeclareTypes();
  dddctrl.allTypesDefined = FALSE;

  DomInitParallel(TypeBndP,TypeBndS);

  ddd_IfInit();

  /* check for correct initialization */
  if ((err=CheckInitParallel())!=0)
  {
    SetHiWrd(err,__LINE__);
    return(err);
  }

  return 0;          /* no error */
}



/****************************************************************************/
/*
   ExitParallel - exit the parallel application on ddd level

   SYNOPSIS:
   int ExitParallel (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function exits the parallel application on ddd level

   RETUR
   N VALUE:
   void
 */
/****************************************************************************/

int ExitParallel (void)
{
  /* free memory allocated by InitParallel */
  if (dddctrl._context!=NULL)
    free(dddctrl._context);

  DDD_Exit();

  return 0;          /* no error */
}

#endif /* ModelP */
