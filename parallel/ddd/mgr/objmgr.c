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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "dddi.h"


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
#define ProcFromId(n)  ((n)&((1<<MAX_PROCBITS_IN_GID)-1))
#define CountFromId(n) (((n)-((n)&((1<<MAX_PROCBITS_IN_GID)-1)))>>MAX_PROCBITS_IN_GID)




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


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

extern int theIdCount;



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

  if (nObjs==0)
    return(NULL);

  locObjs = (DDD_HDR *) AllocTmp(nObjs*sizeof(DDD_HDR));
  if (locObjs==NULL) {
    DDD_PrintError('E', 2210,  STR_NOMEM " in LocalObjectsList");
    return(NULL);
  }

  memcpy(locObjs, ddd_ObjTable, nObjs*sizeof(DDD_HDR));
  qsort(locObjs, nObjs, sizeof(DDD_HDR), sort_ObjListGID);

  return(locObjs);
}


DDD_HDR *LocalCoupledObjectsList (void)
{
  DDD_HDR   *locObjs;

  if (NCPL_GET==0)
    return(NULL);

  locObjs = (DDD_HDR *) AllocTmp(NCPL_GET*sizeof(DDD_HDR));
  if (locObjs==NULL) {
    DDD_PrintError('E', 2211, STR_NOMEM " in LocalCoupledObjectsList");
    return(NULL);
  }

  memcpy(locObjs, ddd_ObjTable, NCPL_GET*sizeof(DDD_HDR));
  qsort(locObjs, NCPL_GET, sizeof(DDD_HDR), sort_ObjListGID);

  return(locObjs);
}


/****************************************************************************/

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

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)

DDD_OBJ DDD_ObjNew (size_t size, DDD_TYPE typ, DDD_PRIO prio, DDD_ATTR attr)
{
  DDD_OBJ obj;

  /* check input parameters */
  if (prio>=MAX_PRIO)
  {
    sprintf(cBuffer,
            "priority must be less than %d in DDD_ObjNew", MAX_PRIO);
    DDD_PrintError('E', 2205, cBuffer);
    HARD_EXIT;
  }
  if (typ>=MAX_TYPEDESC)
  {
    sprintf(cBuffer,
            "DDD-type must be less than %d in DDD_ObjNew", MAX_TYPEDESC);
    DDD_PrintError('E', 2206, cBuffer);
    HARD_EXIT;
  }

  /* get object memory */
  obj = (DDD_OBJ) AllocObj(size, typ, prio, attr);
  if (obj==NULL) {
    DDD_PrintError('E', 2200, STR_NOMEM " in DDD_ObjNew");
    return(NULL);
  }

#       ifdef DebugCreation
  sprintf(cBuffer, "%4d: DDD_ObjNew(size=%d, typ=%d, prio=%d, attr=%d),"
          " ADR=%08x\n",
          me, size, typ, prio, attr, obj);
  DDD_PrintDebug(cBuffer);
#       endif

  return(obj);
}

#else

DDD_OBJ DDD_ObjNew (size_t size, DDD_TYPE typ, DDD_PRIO prio, DDD_ATTR attr)

{
  TYPE_DESC *desc = &(theTypeDefs[typ]);
  DDD_OBJ obj;

  if (desc->handlerALLOCOBJ)
  {
    desc->handlerALLOCOBJ(&obj);
  }
  else
  {
    if (desc->nextFree == desc->arraySize)
    {
      DDD_PrintError('E', 2200, "No more objects in DDD_ObjNew");
      return(0);
    }
    else
      obj = desc->nextFree;

    do
    {
      desc->nextFree++;

      /* cycle around */
      if (desc->nextFree == desc->arraySize) desc->nextFree = 0;

      if (desc->nextFree == obj)
      {
        /* No more free objects */
        desc->nextFree = desc->arraySize;
        break;
      }
    } while (!IsHdrInvalid (&desc->hdr[desc->nextFree]) );
  }

  if (obj != 0)
  {
    OBJ_ATTR(OBJ2HDR(obj,desc)) = attr;
    OBJ_PRIO(OBJ2HDR(obj,desc)) = prio;

#               ifdef DebugCreation
    sprintf(cBuffer, "%d: Allocated obj %d of type %d\n", me, obj, typ);
    DDD_PrintDebug(cBuffer);
    fflush(stdout);
#               endif
  }

  return (obj);
}

#endif

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

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)

void DDD_ObjDelete (DDD_OBJ obj, size_t size, DDD_TYPE typ)
{
  FreeObj((void *)obj, size, typ);
}

#else

void DDD_ObjDelete (DDD_OBJ obj, size_t size, DDD_TYPE typ)

{
  TYPE_DESC *desc = &(theTypeDefs[typ]);

  if (desc->handlerFREEOBJ)
  {
    desc->handlerFREEOBJ(&obj);
  }

  DDD_HdrDestructor (OBJ2HDR(obj,desc));
}

#       endif

/****************************************************************************/
/*                                                                          */
/* Function:  DDD_HdrConstructor                                            */
/*                                                                          */
/****************************************************************************/

/**
        Initiate object's DDD Header.
        This function registers a DDD-object via constructing
        its DDD-header. Each DDD-object is given a unique {\em global ID},
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

   @param hdr   Pointer to the DDD Header which should be constructed
   @param typ   DDD Type of DDD Object
   @param prio  Priority of DDD Object
   @param attr  Attribute of DDD Object
 */

#if defined(C_FRONTEND) || defined(F_FRONTEND)
void DDD_HdrConstructor (DDD_HDR hdr,DDD_TYPE typ,DDD_PRIO prio,DDD_ATTR attr)
{
#endif

#ifdef CPP_FRONTEND
// construct as invalid DDD_Object
DDD_Object::DDD_Object (void)
{
  /* invalidate this DDD_HDR */
  MarkHdrInvalid(&_hdr);
}



// construct as valid DDD_Object
DDD_Object::DDD_Object (DDD_TYPE typ, DDD_PRIO prio, DDD_ATTR attr)
{
  Init(typ, prio, attr);
}


void DDD_Object::Init (DDD_TYPE typ, DDD_PRIO prio, DDD_ATTR attr)
{
  DDD_HDR hdr = &_hdr;

  if (! IsHdrInvalid(hdr))
  {
    sprintf(cBuffer,
            "cannot initialize DDD_Object %08x twice in DDD_Object::Init",
            OBJ_GID(hdr));
    DDD_PrintError('E', 2250, cBuffer);
    HARD_EXIT;
  }
#endif

/* check input parameters */
if (prio>=MAX_PRIO)
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
if (nObjs==MAX_OBJ)
{
  /* TODO update docu */
  /* this is a fatal case. we cant register more objects here */
  DDD_PrintError('F', 2220, "no more objects in DDD_HdrConstructor");
  /* TODO one could try to expand the global tables here. */
  HARD_EXIT;
}

/* insert into theObj array */
ddd_ObjTable[nObjs] = hdr;
OBJ_INDEX(hdr) = nObjs;
nObjs++;
        #else
/* if we dont have WithFullObjectTable, pure local objects without
   copies on other processors aren't registered by DDD. Therefore,
   they don't have a valid OBJ_INDEX field. */
MarkHdrLocal(hdr);
        #endif


/* init object header with defaults */
OBJ_TYPE(hdr)  = typ;
OBJ_PRIO(hdr)  = prio;
OBJ_ATTR(hdr)  = attr;
OBJ_FLAGS(hdr) = 0;

/* create unique GID */
OBJ_GID(hdr)   = MakeUnique(theIdCount++);

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
        "typ=%d, prio=%d, attr=%d), "
        "GID=%08x  INDEX=%d\n",
        me, hdr, typ, prio, attr, OBJ_GID(hdr), OBJ_INDEX(hdr));
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
        via destructing its DDD-header.
        {\em Note:} The DDD-object will be destructed, but its copies
        on remote processors will not be informed by \funk{HdrDestructor}.
        There are two consistent possibilities to delete DDD-objects which
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

#if defined(C_FRONTEND) || defined(F_FRONTEND)
/* in F_FRONTEND, DDD_HdrDestructor is used internally */
void DDD_HdrDestructor (DDD_HDR hdr)
{
#endif
#ifdef CPP_FRONTEND
DDD_Object::~DDD_Object (void)
{
  DDD_HDR hdr = &_hdr;
#endif
COUPLING   *cpl;
int objIndex, xfer_active = XferActive();

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
  XferRegisterDelete(hdr);


objIndex = OBJ_INDEX(hdr);

if (objIndex<NCPL_GET)
{
  /* this is an object with couplings */
  cpl = theCpl[objIndex];

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

  NCPL_DECREMENT;
  nObjs--;

  /* fill slot of deleted obj with last cpl-obj */
  ddd_ObjTable[objIndex] = ddd_ObjTable[NCPL_GET];
  theCpl[objIndex] = theCpl[NCPL_GET];
  theCplN[objIndex] = theCplN[NCPL_GET];
  OBJ_INDEX(ddd_ObjTable[objIndex]) = objIndex;

                #ifdef WithFullObjectTable
  /* fill slot of last cpl-obj with last obj */
  if (NCPL_GET<nObjs)
  {
    ddd_ObjTable[NCPL_GET] = ddd_ObjTable[nObjs];
    OBJ_INDEX(ddd_ObjTable[NCPL_GET]) = NCPL_GET;
  }
                #else
  assert(NCPL_GET==nObjs);
                #endif

  /* dispose all couplings */
  DisposeCouplingList(cpl);
}
else
{
                #ifdef WithFullObjectTable
  /* this is an object without couplings */
  /* deletion is not dangerous (no consistency problem) */
  nObjs--;

  /* fill slot of deleted obj with last obj */
  ddd_ObjTable[objIndex] = ddd_ObjTable[nObjs];
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

#ifdef F_FRONTEND
void DDD_ObjGet (DDD_TYPE *ftyp, DDD_PRIO *fprio, DDD_ATTR *fattr,
                 DDD_OBJ *fobj)
{
  TYPE_DESC *desc = &(theTypeDefs[*ftyp]);

  /* get the index of the object */
  *fobj = (DDD_OBJ) DDD_ObjNew (0, *ftyp, *fprio, *fattr);

  if (*fobj==0)
  {
    DDD_PrintError('E', 2200, STR_NOMEM " in DDD_ObjGet");
    return;
  }

  /* call DDD_HdrConstructor */
  DDD_HdrConstructor (OBJ2HDR(*fobj,desc), *ftyp, *fprio, *fattr);

  return;
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

#if defined(F_FRONTEND)
void DDD_ObjUnGet (DDD_OBJ *fobj, DDD_TYPE *ftyp)

{
  DDD_OBJ obj  = *fobj;
  TYPE_DESC *desc = &(theTypeDefs[*ftyp]);
  DDD_HDR hdr  = OBJ2HDR (obj, desc);

  /* call DDD_HDR-destructor */
  DDD_HdrDestructor (hdr);

  /* free raw memory */
  DDD_ObjDelete (obj, 0, *ftyp);
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
  if (nObjs==MAX_OBJ)
  {
    /* TODO update docu */
    /* this is a fatal case. we cant register more objects here */
    DDD_PrintError('F', 2220, "no more objects in DDD_HdrConstructorCopy");
    /* TODO one could try to expand the global tables here. */
  }

  /* insert into theObj array */
  ddd_ObjTable[nObjs] = newhdr;
  OBJ_INDEX(newhdr) = nObjs;
  nObjs++;
        #else
  MarkHdrLocal(newhdr);
  assert(nObjs==NCPL_GET);
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
  if (objIndex<NCPL_GET)
    ddd_ObjTable[objIndex] = newhdr;
        #endif

  /* change pointers from couplings to object */
  if (objIndex<NCPL_GET)
  {
    COUPLING *cpl = theCpl[objIndex];

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

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)

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

#else

static void MsgToObj (TYPE_DESC *desc, DDD_OBJ target, char *source)
{
  DDD_HDR hdr;
  int i, i2;

  /* Copy the header */
  hdr = &(desc->hdr[target]);
  hdr->typ   = (unsigned char) *source++;
  hdr->prio  = (unsigned char) *source++;
  hdr->attr  = (unsigned char) *source++;
  hdr->flags = (unsigned char) *source++;
  source    += sizeof (unsigned int);           /* myIndex is local */
  hdr->gid   = *((unsigned int *)source); source += sizeof (unsigned int);
  source    += 4;                                               /* is local */

  /* Copy the data */
  for (i = 0; i < desc->nElements; i++)
  {
    ELEM_DESC *elem = &(desc->element[i]);
    char      *d    = elem->array + elem->size * target;

    if (elem->type != EL_LDATA)
      for (i2 = elem->size; i2; i2--) *d++ = *source++;
  }
}


void ObjCopyGlobalData (TYPE_DESC *desc, DDD_OBJ target, DDD_OBJ source,
                        size_t size)
{
  MsgToObj (desc, target, (char *) source);

#       ifdef DebugCreation
  sprintf(cBuffer, "%4d: ObjCopyGlobalData(%08x <- %08x),"
          " TYP=%d  GID=%08x  INDEX=%d\n",
          me, OBJ2HDR(target,desc), source,
          OBJ_TYPE(OBJ2HDR(target,desc)),
          OBJ_GID(OBJ2HDR(target,desc)),
          OBJ_INDEX(OBJ2HDR(target,desc)));
  DDD_PrintDebug(cBuffer);
#       endif
}

#endif


/****************************************************************************/

#ifdef C_FRONTEND
DDD_HDR DDD_SearchHdr (DDD_GID gid)
{
#endif
#ifdef CPP_FRONTEND
DDD_HDR DDD_Library::SearchHdr (DDD_GID gid)
{
#endif
#ifdef F_FRONTEND
DDD_HDR DDD_SearchHdr (DDD_GID *_gid)
{
  DDD_GID gid = *_gid;
#endif
int i;

i=0;
while (i<nObjs && OBJ_GID(ddd_ObjTable[i])!=gid)
  i++;

if (i<nObjs)
{
  return(ddd_ObjTable[i]);
}
else
  return(NULL);
}


/****************************************************************************/
