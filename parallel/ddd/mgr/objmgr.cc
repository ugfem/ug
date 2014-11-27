// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      objmgr.c                                                      */
/*                                                                          */
/* Purpose:   creation and deletion of objects                              */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/02/21 kb  begin                                            */
/*            95/11/03 kb  complete redesign of objmgr-interface, C++ style */
/*            95/11/15 kb  arbitrary offset for DDD_HEADER                  */
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
#include <config.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>

#include "dddi.h"

USING_UG_NAMESPACES

/* PPIF namespace: */
USING_PPIF_NAMESPACE

START_UGDIM_NAMESPACE

/*
   #define DebugCreation
   #define DebugDeletion
 */


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/


#define MakeUnique(n)  (((n)<<MAX_PROCBITS_IN_GID)+me)
#define ProcFromId(n)  ((n)& ((1<<MAX_PROCBITS_IN_GID)-1))
#define CountFromId(n) (((n)-((n)& ((1<<MAX_PROCBITS_IN_GID)-1)))>>MAX_PROCBITS_IN_GID)




/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)



/* local unique ID count */
static DDD_GID theIdCount;



/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/




/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/



static int sort_ObjListGID (const void *e1, const void *e2)
{
  DDD_HDR o1, o2;

  o1 = *((DDD_HDR *)e1);
  o2 = *((DDD_HDR *)e2);

  if (OBJ_GID(o1) < OBJ_GID(o2)) return(-1);
  if (OBJ_GID(o1) == OBJ_GID(o2)) return(0);
  return(1);
}


DDD_HDR *LocalObjectsList (void)
{
  DDD_HDR   *locObjs;

  if (ddd_nObjs==0)
    return(NULL);

  locObjs = (DDD_HDR *) AllocTmpReq (ddd_nObjs*sizeof(DDD_HDR), TMEM_OBJLIST);
  if (locObjs==NULL) {
    DDD_PrintError('E', 2210,  STR_NOMEM " in LocalObjectsList");
    return(NULL);
  }

  memcpy(locObjs, ddd_ObjTable, ddd_nObjs*sizeof(DDD_HDR));
  qsort(locObjs, ddd_nObjs, sizeof(DDD_HDR), sort_ObjListGID);

  return(locObjs);
}

void FreeLocalObjectsList (DDD_HDR *locObjs)
{
  if (locObjs==NULL)
    return;

  FreeTmpReq(locObjs, ddd_nObjs*sizeof(DDD_HDR), TMEM_OBJLIST);
}



DDD_HDR *LocalCoupledObjectsList (void)
{
  DDD_HDR   *locObjs;

  if (NCpl_Get==0)
    return(NULL);

  locObjs = (DDD_HDR *) AllocTmpReq (NCpl_Get*sizeof(DDD_HDR), TMEM_OBJLIST);
  if (locObjs==NULL) {
    DDD_PrintError('E', 2211, STR_NOMEM " in LocalCoupledObjectsList");
    return(NULL);
  }

  memcpy(locObjs, ddd_ObjTable, NCpl_Get*sizeof(DDD_HDR));
  qsort(locObjs, NCpl_Get, sizeof(DDD_HDR), sort_ObjListGID);

  return(locObjs);
}


void FreeLocalCoupledObjectsList (DDD_HDR *locObjs)
{
  if (locObjs==NULL)
    return;

  FreeTmpReq(locObjs, NCpl_Get*sizeof(DDD_HDR), TMEM_OBJLIST);
}



/****************************************************************************/


void ddd_EnsureObjTabSize (int n)
{
  DDD_HDR *old_ObjTable   = ddd_ObjTable;
  int old_ObjTabSize = ddd_ObjTabSize;

  /* if size is large enough, we are already finished. */
  if (old_ObjTabSize >= n)
    return;

  /* set new size */
  ddd_ObjTabSize = n;

  /* allocate new object table */
  ddd_ObjTable = (DDD_HDR *) AllocTmp(sizeof(DDD_HDR) * ddd_ObjTabSize);
  if (ddd_ObjTable==NULL)
  {
    sprintf(cBuffer, STR_NOMEM " for object table of size %ld",
            ((long)ddd_ObjTabSize) * sizeof(DDD_HDR));
    DDD_PrintError('E', 2223, cBuffer);
    HARD_EXIT;
  }

  /* copy data from old cpl-table to new one, assuming the old one is full */
  memcpy(ddd_ObjTable, old_ObjTable, sizeof(DDD_HDR) * old_ObjTabSize);

  /* free old one */
  FreeTmp(old_ObjTable,0);

  /* issue a warning in order to inform user */
  sprintf(cBuffer, "increased object table, now %d entries", ddd_ObjTabSize);
  DDD_PrintError('W', 2224, cBuffer);
}


/****************************************************************************/

/*
        Description of ObjMgr Interfaces

        Raw-Memory-Interface:    DDD_ObjNew, DDD_ObjDelete
        Constructor-Interface:   DDD_HdrConstructor, DDD_HdrDestructor
        Application-Interface:   DDD_Get, DDD_UnGet
 */


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_ObjNew                                                    */
/*                                                                          */
/****************************************************************************/

/* Purpose:   get raw memory for new DDD-object.                            */
/*                                                                          */
/* Input:     size:  memory size of object                                  */
/*            typ:   DDD_TYPE of object                                     */
/*            prio:  DDD_PRIO of object                                     */
/*            attr:  attribute of distributed object                        */
/*                                                                          */
/* Output:    pointer to raw memory                                         */
/*                                                                          */
/****************************************************************************/

/**
        Allocate raw memory for new \ddd{object}.
        This function dynamically creates raw memory for a new \ddd{object}.
        Therefore, the user-supplied memory manager function \memmgrfunk{AllocOMEM}
        is called to allocate the necessary memory. Although the caller must
        supply the object's priority and attribute, its header will not be
        initialized by \funk{ObjNew}; the parameters are used for smart
        memory allocation, only.

        The function \funk{ObjNew} and \funk{ObjDelete}, its corresponding
        deletion function, form the ObjManager's {\em raw memory interface}.

        DDD users who use the #C_FRONTEND# may prefer the more
        elaborate {\em application interface}, consisting of the functions
        \funk{ObjGet} and \funk{ObjUnGet}. DDD users who use the language
        C++ (object-oriented style) with #CPP_FRONTEND# will use the
        {\em raw memory interface} together with the
        {\em constructor/destructor interface} (\funk{HdrConstructor},
        \funk{HdrDestructor}, \funk{HdrConstructorMove})
        in order to integrate DDD into C++ style object management easily.

        For variable-sized \ddd{objects}, the parameter {\em aSize}
        may differ from the size specified during the corresponding
        \funk{TypeDefine}-call.

   @return pointer to free memory block for the \ddd{object}
   @param  aSize   memory size of the new object
   @param  aType   \ddd{type} of the new object
   @param  aPrio   \ddd{priority} of the new object
   @param  aAttr   \ddd{attribute} of the new object
 */


DDD_OBJ DDD_ObjNew (size_t aSize, DDD_TYPE aType,
                    DDD_PRIO aPrio, DDD_ATTR aAttr)
{
  DDD_OBJ obj;

  /* check input parameters */
  if (aPrio>=MAX_PRIO)
  {
    sprintf(cBuffer,
            "priority must be less than %d in DDD_ObjNew", MAX_PRIO);
    DDD_PrintError('E', 2205, cBuffer);
    HARD_EXIT;
  }
  if (aType>=MAX_TYPEDESC)
  {
    sprintf(cBuffer,
            "DDD-type must be less than %d in DDD_ObjNew", MAX_TYPEDESC);
    DDD_PrintError('E', 2206, cBuffer);
    HARD_EXIT;
  }

  /* get object memory */
  obj = (DDD_OBJ) AllocObj(aSize, aType, aPrio, aAttr);
  if (obj==NULL) {
    DDD_PrintError('E', 2200, STR_NOMEM " in DDD_ObjNew");
    return(NULL);
  }

#       ifdef DebugCreation
  sprintf(cBuffer, "%4d: DDD_ObjNew(aSize=%d, type=%d, prio=%d, attr=%d),"
          " ADR=%08x\n",
          me, aSize, aType, aPrio, aAttr, obj);
  DDD_PrintDebug(cBuffer);
#       endif

  return(obj);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_ObjDelete                                                 */
/*                                                                          */
/* Purpose:   free raw memory from DDD-object                               */
/*                                                                          */
/* Input:     obj:   object header address                                  */
/*            size:  memory size of object                                  */
/*            typ:   DDD_TYPE of object                                     */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/


void DDD_ObjDelete (DDD_OBJ obj, size_t size, DDD_TYPE typ)
{
  FreeObj((void *)obj, size, typ);
}


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_HdrConstructor                                            */
/*                                                                          */
/****************************************************************************/

/**
        Initiate object's \ddd{header}.
        This function registers a \ddd{object} via constructing
        its DDD-header. Each \ddd{object} is given a unique {\em global ID},
        which is stored in the DDD-header together with the object's
        properties (type\_id, prio, attr) and additional data used by DDD.

        The function \funk{HdrConstructor} and its corresponding deletion
        function \funk{HdrDestructor} form the ObjManager's
        {\em constructor/destructor interface}. This interface will be
        useful for C++ users, together with the {\em raw memory interface}
        functions \funk{ObjNew} and \funk{ObjDelete}.
        DDD users who use the language C may prefer the more
        elaborate {\em application interface}, consisting of the functions
        \funk{ObjGet} and \funk{ObjUnGet}.

        Due to the construction of global IDs, an overflow error may occur
        after many calls to \funk{HdrConstructor}
        (\MaxUniqueGids\ times in \Version).

   @param aHdr   pointer to the \ddd{header} which should be constructed
   @param aType  \ddd{type} of \ddd{object}
   @param aPrio  \ddd{priority} of \ddd{object}
   @param aAttr  \ddd{attribute} of \ddd{object}
 */

#if defined(C_FRONTEND)
void DDD_HdrConstructor (DDD_HDR aHdr, DDD_TYPE aType,
                         DDD_PRIO aPrio, DDD_ATTR aAttr)
{
#endif

#ifdef CPP_FRONTEND
// construct as invalid DDD_Object
DDD_Object::DDD_Object (void)
{
  /* invalidate this DDD_HDR */
  MarkHdrInvalid(this);
}



// construct as valid DDD_Object
DDD_Object::DDD_Object (DDD_TYPE aType, DDD_PRIO aPrio, DDD_ATTR aAttr)
{
  Init(aType, aPrio, aAttr);
}


void DDD_Object::Init (DDD_TYPE aType, DDD_PRIO aPrio, DDD_ATTR aAttr)
{
  DDD_HDR aHdr = this;

  if (! IsHdrInvalid(aHdr))
  {
    sprintf(cBuffer,
            "cannot initialize DDD_Object %08x twice in DDD_Object::Init",
            OBJ_GID(aHdr));
    DDD_PrintError('E', 2250, cBuffer);
    HARD_EXIT;
  }
#endif

/* check input parameters */
if (aPrio>=MAX_PRIO)
{
  sprintf(cBuffer,
          "priority must be less than %d in DDD_HdrConstructor", MAX_PRIO);
  DDD_PrintError('E', 2225, cBuffer);
  HARD_EXIT;
}

        #ifdef WithFullObjectTable
/* in case of FullObjectTable, we register each header in the
   global ddd_ObjTable. */

/* check whether there are available objects */
if (ddd_nObjs==ddd_ObjTabSize)
{
  /* TODO update docu */
  /* this is a fatal case. we cant register more objects here */
  DDD_PrintError('F', 2220, "no more objects in DDD_HdrConstructor");
  /* TODO one could try to expand the global tables here. */
  HARD_EXIT;
}

/* insert into theObj array */
ddd_ObjTable[ddd_nObjs] = aHdr;
OBJ_INDEX(aHdr) = ddd_nObjs;
ddd_nObjs++;
        #else
/* if we dont have WithFullObjectTable, pure local objects without
   copies on other processors aren't registered by DDD. Therefore,
   they don't have a valid OBJ_INDEX field. */
MarkHdrLocal(aHdr);
        #endif


/* init object header with defaults */
OBJ_TYPE(aHdr)  = aType;
OBJ_PRIO(aHdr)  = aPrio;
OBJ_ATTR(aHdr)  = aAttr;
OBJ_FLAGS(aHdr) = 0;

/* create unique GID */
OBJ_GID(aHdr)   = MakeUnique(theIdCount++);

/* check overflow of global id numbering */
if (MakeUnique(theIdCount) <= MakeUnique(theIdCount-1))
{
  /* TODO update docu */
  DDD_PrintError('F', 2221, "global ID overflow DDD_HdrConstructor");
  /* TODO one could try to renumber all objects here. */
  HARD_EXIT;
}

#       ifdef DebugCreation
sprintf(cBuffer, "%4d: DDD_HdrConstructor(adr=%08x, "
        "type=%d, prio=%d, attr=%d), "
        "GID=%08x  INDEX=%d\n",
        me, aHdr, aType, aPrio, aAttr, OBJ_GID(aHdr), OBJ_INDEX(aHdr));
DDD_PrintDebug(cBuffer);
#       endif
}



#ifdef CPP_FRONTEND
DDD_IndexObject::DDD_IndexObject (DDD_TYPE typ,
                                  DDD_INDEX idx, DDD_PRIO prio, DDD_ATTR attr)
  : DDD_Object(typ, prio, attr)
{
  _index = idx;
}
#endif



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_HdrDestructor                                             */
/*                                                                          */
/****************************************************************************/

/**
        Remove object's header from DDD management.
        This function removes an object from DDD-management
        via destructing its \ddd{header}.
        {\em Note:} The \ddd{object} will be destroyed, but its copies
        on remote processors will not be informed by \funk{HdrDestructor}.
        There are two consistent possibilities to delete \ddd{objects} which
        have copies on remote processors:

        \begin{itemize}
        \item In order to delete only this local object copy, use
        \funk{XferDelete} during a DDD Transfer-operation. This will
        inform the remote copies of the deletion.
        \item In order to delete a distributed object (\ie, all its
        object copies), use function \funk{ObjUnGet} (when using the
        {\em application interface})
        or a combination of functions
        \funk{HdrDestructor} / \funk{ObjDelete} (when not using
        the {\em application interface}) for all copies.
        \end{itemize}

        The function \funk{HdrDestructor} and its corresponding creation
        function \funk{HdrConstructor} form the ObjManager's
        {\em constructor/destructor interface}. This interface will be
        useful for C++ users, together with the {\em raw memory interface}
        functions \funk{ObjNew} and \funk{ObjDelete}.
        DDD users who use the language C may prefer the more
        elaborate {\em application interface}, consisting of the functions
        \funk{ObjGet} and \funk{ObjUnGet}.

   @param hdr  the object's DDD Header
 */

#if defined(C_FRONTEND)
void DDD_HdrDestructor (DDD_HDR hdr)
{
#endif
#ifdef CPP_FRONTEND
DDD_Object::~DDD_Object (void)
{
  DDD_HDR hdr = this;
#endif
COUPLING   *cpl;
int objIndex, xfer_active = ddd_XferActive();

#       ifdef DebugDeletion
sprintf(cBuffer, "%4d: DDD_HdrDestructor(adr=%08x, "
        "typ=%d, prio=%d, attr=%d), "
        "GID=%08x  INDEX=%d\n",
        me, hdr, OBJ_TYPE(hdr), OBJ_PRIO(hdr), OBJ_ATTR(hdr),
        OBJ_GID(hdr), OBJ_INDEX(hdr));
DDD_PrintDebug(cBuffer);
#       endif


if (IsHdrInvalid(hdr))
{
  /* DDD_HDR is invalid, so destructor is useless */
  return;
}

/* formally, the object's GID should be returned here */


/* if currently in xfer, register deletion for other processors */
if (xfer_active)
  ddd_XferRegisterDelete(hdr);


objIndex = OBJ_INDEX(hdr);

if (objIndex<NCpl_Get)
{
  /* this is an object with couplings */
  cpl = IdxCplList(objIndex);

  /* if not during xfer, deletion may be inconsistent */
  if (!xfer_active)
  {
    /* deletion is dangerous, distributed object might get
       inconsistent. */
    if (DDD_GetOption(OPT_WARNING_DESTRUCT_HDR)==OPT_ON)
    {
      sprintf(cBuffer,
              "inconsistency by deleting gid=%08x in DDD_HdrDestructor",
              OBJ_GID(hdr));
      DDD_PrintError('W', 2230, cBuffer);
    }
  }

  NCpl_Decrement;
  ddd_nObjs--;

  /* fill slot of deleted obj with last cpl-obj */
  ddd_ObjTable[objIndex] = ddd_ObjTable[NCpl_Get];
  IdxCplList(objIndex) = IdxCplList(NCpl_Get);
  IdxNCpl(objIndex) = IdxNCpl(NCpl_Get);
  OBJ_INDEX(ddd_ObjTable[objIndex]) = objIndex;

                #ifdef WithFullObjectTable
  /* fill slot of last cpl-obj with last obj */
  if (NCpl_Get<ddd_nObjs)
  {
    ddd_ObjTable[NCpl_Get] = ddd_ObjTable[ddd_nObjs];
    OBJ_INDEX(ddd_ObjTable[NCpl_Get]) = NCpl_Get;
  }
                #else
  assert(NCpl_Get==ddd_nObjs);
                #endif

  /* dispose all couplings */
  DisposeCouplingList(cpl);
}
else
{
                #ifdef WithFullObjectTable
  /* this is an object without couplings */
  /* deletion is not dangerous (no consistency problem) */
  ddd_nObjs--;

  /* fill slot of deleted obj with last obj */
  ddd_ObjTable[objIndex] = ddd_ObjTable[ddd_nObjs];
  OBJ_INDEX(ddd_ObjTable[objIndex]) = objIndex;
                #endif
}

/* invalidate this DDD_HDR */
MarkHdrInvalid(hdr);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_ObjGet                                                    */
/*                                                                          */
/* Purpose:   get new DDD object for a given DDD_TYPE                       */
/*                                                                          */
/*            DDD_ObjGet allows different size for each object instance,    */
/*            differing from the size computed during DDD_TypeDefine        */
/*                                                                          */
/* Input:     typ: DDD_TYPE of object                                       */
/*                                                                          */
/* Output:    pointer to memory, which is raw except from the               */
/*              constructed DDD_HDR.                                        */
/*                                                                          */
/****************************************************************************/

#if defined(C_FRONTEND)
DDD_OBJ DDD_ObjGet (size_t size, DDD_TYPE typ, DDD_PRIO prio, DDD_ATTR attr)
{
  DDD_OBJ obj;
  TYPE_DESC  *desc = &(theTypeDefs[typ]);

  /* check input parameters */
  if (prio<0 || prio>=MAX_PRIO)
  {
    sprintf(cBuffer,
            "priority must be less than %d in DDD_ObjGet", MAX_PRIO);
    DDD_PrintError('E', 2235, cBuffer);
    HARD_EXIT;
  }

  /* get raw memory */
  obj = (DDD_OBJ) DDD_ObjNew(size, typ, prio, attr);
  if (obj==NULL) {
    DDD_PrintError('E', 2200, STR_NOMEM " in DDD_ObjGet");
    return(NULL);
  }

  if ((desc->size != size) && (DDD_GetOption(OPT_WARNING_VARSIZE_OBJ)==OPT_ON))
  {
    DDD_PrintError('W', 2200,
                   "object size differs from declared size in DDD_ObjGet");
  }

  if ((desc->size > size) && (DDD_GetOption(OPT_WARNING_SMALLSIZE)==OPT_ON))
  {
    DDD_PrintError('W', 2201,
                   "object size smaller than declared size in DDD_ObjGet");
  }


  /* call DDD_HdrConstructor */
  DDD_HdrConstructor(OBJ2HDR(obj,desc), typ, prio, attr);

  return(obj);
}
#endif



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_ObjUnGet                                                  */
/*                                                                          */
/* Purpose:   remove object from DDD management and free its memory         */
/*                                                                          */
/* Input:     obj: object header address                                    */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#if defined(C_FRONTEND)
void DDD_ObjUnGet (DDD_HDR hdr, size_t size)

{
  DDD_TYPE typ = OBJ_TYPE(hdr);
  TYPE_DESC  *desc = &(theTypeDefs[typ]);
  DDD_OBJ obj = HDR2OBJ(hdr,desc);

  if ((desc->size != size) && (DDD_GetOption(OPT_WARNING_VARSIZE_OBJ)==OPT_ON))
  {
    DDD_PrintError('W', 2299,
                   "object size differs from declared size in DDD_ObjUnGet");
  }

  /* call DDD_HDR-destructor */
  DDD_HdrDestructor(hdr);

  /* free raw memory */
  DDD_ObjDelete(obj, size, typ);
}
#endif



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_HdrConstructorCopy                                        */
/*                                                                          */
/* Purpose:   create DDD_HDR copy from message original                     */
/*                                                                          */
/* Input:     newhdr: new DDD_HDR                                           */
/*            prio:   DDD_PRIO of new copy                                  */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void DDD_HdrConstructorCopy (DDD_HDR newhdr, DDD_PRIO prio)
{
  /* check input parameters */
  if (prio>=MAX_PRIO)
  {
    sprintf(cBuffer,
            "priority must be less than %d in DDD_HdrConstructorCopy", MAX_PRIO);
    DDD_PrintError('E', 2245, cBuffer);
    HARD_EXIT;
  }

        #ifdef WithFullObjectTable
  /* check whether there are available objects */
  if (ddd_nObjs==ddd_ObjTabSize)
  {
    /* TODO update docu */
    /* this is a fatal case. we cant register more objects here */
    DDD_PrintError('F', 2220, "no more objects in DDD_HdrConstructorCopy");
    /* TODO one could try to expand the global tables here. */
  }

  /* insert into theObj array */
  ddd_ObjTable[ddd_nObjs] = newhdr;
  OBJ_INDEX(newhdr) = ddd_nObjs;
  ddd_nObjs++;
        #else
  MarkHdrLocal(newhdr);
  assert(ddd_nObjs==NCpl_Get);
        #endif

  /* init LDATA components. GDATA components will be copied elsewhere */
  OBJ_PRIO(newhdr)  = prio;


#       ifdef DebugCreation
  sprintf(cBuffer, "%4d: DDD_HdrConstructorCopy(adr=%08x, prio=%d), "
          "GID=%08x  INDEX=%d  ATTR=%d\n",
          me, newhdr, prio, OBJ_GID(newhdr), OBJ_INDEX(newhdr), OBJ_ATTR(newhdr));
  DDD_PrintDebug(cBuffer);
#       endif
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_HdrConstructorMove                                        */
/*                                                                          */
/* Purpose:   create DDD_HDR copy inside local memory,                      */
/*            simultaneously destruct original DDD_HDR                      */
/*                                                                          */
/* Input:     newhdr: new DDD_HDR                                           */
/*            oldhdr: old DDD_HDR                                           */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void DDD_HdrConstructorMove (DDD_HDR newhdr, DDD_HDR oldhdr)
{
  int objIndex = OBJ_INDEX(oldhdr);


  /* copy all components */
  OBJ_INDEX(newhdr) = OBJ_INDEX(oldhdr);
  OBJ_TYPE(newhdr)  = OBJ_TYPE(oldhdr);
  OBJ_PRIO(newhdr)  = OBJ_PRIO(oldhdr);
  OBJ_ATTR(newhdr)  = OBJ_ATTR(oldhdr);
  OBJ_FLAGS(newhdr) = OBJ_FLAGS(oldhdr);
  OBJ_GID(newhdr)   = OBJ_GID(oldhdr);


  /* change all references from DDD to oldhdr */

  /* change entry of theObj array */
        #ifdef WithFullObjectTable
  ddd_ObjTable[objIndex] = newhdr;
        #else
  if (objIndex<NCpl_Get)
    ddd_ObjTable[objIndex] = newhdr;
        #endif

  /* change pointers from couplings to object */
  if (objIndex<NCpl_Get)
  {
    COUPLING *cpl = IdxCplList(objIndex);

    for(; cpl!=NULL; cpl=CPL_NEXT(cpl)) {
      cpl->obj = newhdr;
    }

    /* invalidate update obj-shortcut tables from IF module */
    IFInvalidateShortcuts(OBJ_TYPE(newhdr));
  }

  /* invalidate old DDD_HDR */
  MarkHdrInvalid(oldhdr);
}



/****************************************************************************/
/*                                                                          */
/* Function:  ObjCopyGlobalData                                             */
/*                                                                          */
/* Purpose:   copy DDD-object from message to memory. elements marked       */
/*            EL_LDATA (local data) will not be copied. this is done        */
/*            in an efficient way by using the DDD_TYPEs copy-mask (which   */
/*            has been set up during StructRegister()).                     */
/*                                                                          */
/*            CopyByMask() is a support function doing the actual work.     */
/*            This function could be more efficient by copying more than    */
/*            one byte at a time.                                           */
/*                                                                          */
/* Input:     target: DDD_OBJ address of target memory                      */
/*            source: DDD_OBJ address of source object                      */
/*            size:   size of object in bytes (might be != desc->size)      */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/


static void CopyByMask (TYPE_DESC *desc, DDD_OBJ target, DDD_OBJ source)
{
  unsigned char *maskp;
  unsigned char *s=(unsigned char *)source, *t=(unsigned char *)target;
  int i;

#       ifdef DebugCreation
  sprintf(cBuffer, "%4d: CopyByMask(%s, size=%d, to=%08x, from=%08x)\n",
          me, desc->name, desc->size, target, source);
  DDD_PrintDebug(cBuffer);
#       endif

  maskp = desc->cmask;

  /* copy all bits set in cmask from source to target */
  for(i=0; i<desc->size; i++)
  {
    unsigned char negmask = *maskp^0xff;
    *t = (*s & *maskp) | (*t & negmask);

    t++; s++; maskp++;              /* Paragon didn't manage postfix increment */
  }
}


void ObjCopyGlobalData (TYPE_DESC *desc,
                        DDD_OBJ target, DDD_OBJ source, size_t size)
{
  /*
          normally size will be equal to desc->size (for fixed-sized objects).
          for variable-sized objects size depends on what the sender put into
          the message
   */

  CopyByMask(desc, target, source);

  /* copy remainder as EL_GDATA */
  if (size > desc->size)
  {
    memcpy(((char *)target)+desc->size,
           ((char *)source)+desc->size, size-desc->size);
  }

#       ifdef DebugCreation
  sprintf(cBuffer, "%4d: ObjCopyGlobalData(%08x <- %08x, size=%d),"
          " TYP=%d  GID=%08x  INDEX=%d\n",
          me, OBJ2HDR(target,desc), source, size,
          OBJ_TYPE(OBJ2HDR(target,desc)),
          OBJ_GID(OBJ2HDR(target,desc)),
          OBJ_INDEX(OBJ2HDR(target,desc)));
  DDD_PrintDebug(cBuffer);
#       endif
}



/****************************************************************************/

#ifdef C_FRONTEND
DDD_HDR DDD_SearchHdr (DDD_GID gid)
{
#endif
#ifdef CPP_FRONTEND
DDD_HDR DDD_Library::SearchHdr (DDD_GID gid)
{
#endif
int i;

i=0;
while (i<ddd_nObjs && OBJ_GID(ddd_ObjTable[i])!=gid)
  i++;

if (i<ddd_nObjs)
{
  return(ddd_ObjTable[i]);
}
else
  return(NULL);
}


/****************************************************************************/


void ddd_ObjMgrInit (void)
{
  /* sanity check: does the DDD_PROC type have enough bits? */
  if (sizeof(DDD_PROC)*8 < MAX_PROCBITS_IN_GID)
    DDD_PrintError('F', 666, "DDD_PROC isn't large enough for MAX_PROCBITS_IN_GID bits");

  theIdCount = 1;        /* start with 1, for debugging reasons */

  /* allocate first (smallest) object table */
  ddd_ObjTable = (DDD_HDR *) AllocTmp(sizeof(DDD_HDR) * MAX_OBJ_START);
  if (ddd_ObjTable==NULL)
  {
    DDD_PrintError('E', 2222, STR_NOMEM " for initial object table");
    HARD_EXIT;
  }

  ddd_ObjTabSize = MAX_OBJ_START;
}


void ddd_ObjMgrExit (void)
{
  FreeTmp(ddd_ObjTable,0);
}


/****************************************************************************/

END_UGDIM_NAMESPACE
