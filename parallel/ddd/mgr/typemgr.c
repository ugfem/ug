// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      typemgr.c                                                     */
/*                                                                          */
/* Purpose:   declaring and defining DDD_TYPEs                              */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/22 kb  begin                                            */
/*            94/09/13 kb  added DDD_Status                                 */
/*            95/11/03 kb  complete rewrite of StructRegister code          */
/*            95/11/16 kb  copied from main.c, introduced clean type concept*/
/*            95/12/07 jb  added fortran frontend                           */
/*            96/04/02 kb  added EL_CONTINUE feature to TypeDefine()        */
/*            97/02/12 kb  added CPP_FRONTEND                               */
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
#include <stdarg.h>
#include <string.h>

#include "dddi.h"


/*
   #define DebugTypeDefine
   #define DebugCopyMask
   #define DebugNoStructCompress
 */


/****************************************************************************/
/*                                                                          */
/* definition of constants                                                  */
/*                                                                          */
/****************************************************************************/

enum DDD_TypeModes
{
  DDD_TYPE_INVALID = 0,            /* DDD_TYPE not declared, not defined      */
  DDD_TYPE_DECLARED,               /* DDD_TYPE declared, but not defined      */
  DDD_TYPE_CONTDEF,                /* DDD_TYPE declared and partially defined */
  DDD_TYPE_DEFINED                 /* DDD_TYPE declared and defined           */
};


#ifdef CPP_FRONTEND
enum DDD_TypeStorageModes
{
  STORAGE_STRUCT,                  /* objects of this type are structs/classes*/
  STORAGE_ARRAY                    /* objects of this type are collections of
                                      contiguous arrays                       */
};
#endif



/* macros for easier switching of FRONTENDs */

#ifdef F_FRONTEND
#define FTYPE  *
#else
#define FTYPE
#endif

#ifdef CPP_FRONTEND
#define CPP_STRUCT(d)     ((d)->storage==STORAGE_STRUCT)
#define CPP_ARRAY(d)      ((d)->storage==STORAGE_ARRAY)
#define CPP_AND           &&
#else
#define CPP_STRUCT(d)
#define CPP_ARRAY(d)
#define CPP_AND
#endif


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/* global table of DDD_TYPE definitions */
TYPE_DESC theTypeDefs[MAX_TYPEDESC];


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only                  */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)




/* overall number of DDD_TYPE definitions */
static int nDescr;



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/

static void InitHandlers (TYPE_DESC *);


int ddd_TypeDefined (TYPE_DESC *desc)
{
  return(desc->mode==DDD_TYPE_DEFINED);
}


/****************************************************************************/


/*
        sort pointers to ELEM_DESC according to their offset
 */

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
static int sort_el_offset (const void *i1, const void *i2)
{
  ELEM_DESC *e1 = (ELEM_DESC *)i1;
  ELEM_DESC *e2 = (ELEM_DESC *)i2;

  if (e1->offset < e2->offset) return(-1);
  if (e1->offset > e2->offset) return(1);
  return(0);
}
#endif


/*
        print out error message during TypeDefine process

        error occurred during TypeDefine for desc, with argument argno
 */

static char *RegisterError (TYPE_DESC *desc, int argno, char *txt)
{
  if (argno==0)
  {
    sprintf(cBuffer, "%s in DDD_TypeDefine(\"%s/%d\")",
            txt, desc->name, desc->currTypeDefCall);
  }
  else
  {
    sprintf(cBuffer, "%s, arg %d of DDD_TypeDefine(\"%s/%d\")",
            txt, argno, desc->name, desc->currTypeDefCall);
  }

  return cBuffer;
}



/*
        check ELEM_DESC for plausibility
 */

static int CheckBounds (TYPE_DESC *desc, ELEM_DESC *el, int argno)
{
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  if (CPP_STRUCT(desc) CPP_AND (el->offset<0))
  {
    DDD_PrintError('E', 9900,
                   RegisterError(desc,argno, "negative offset"));
    return(ERROR);
  }
#endif
#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
  if (CPP_ARRAY(desc) CPP_AND (!el->array))
  {
    DDD_PrintError ('E', 9999,
                    RegisterError(desc,argno, "no array supplied"));
    return(ERROR);
  }
#endif

  if (el->size<=0)
  {
    DDD_PrintError('E', 9901,
                   RegisterError(desc,argno, "illegal element size"));
    return (ERROR);
  }

  return 0;
}



/*
        check ELEM_DESC list of given TYPE_DESC for bad overlapping
 */

static int CheckOverlapEls (TYPE_DESC *desc)
{
  char buf[64];
  int i;
  int ok = TRUE;

  for(i=0; i<desc->nElements; i++)
  {
    ELEM_DESC *e1 = &desc->element[i];

    if (i<desc->nElements-1)
    {
      ELEM_DESC *e2 = &desc->element[i+1];
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      if (CPP_STRUCT(desc) CPP_AND (e1->offset+e1->size > e2->offset))
      {
        ok = FALSE;
        sprintf(buf, "element too big (offset=%d)", e1->offset);
        DDD_PrintError('E', 9902, RegisterError(desc, 0, buf));
      }
#endif
#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
      if (CPP_ARRAY(desc) CPP_AND
            (e1->array+(e1->size*desc->arraySize) > e2->array))
      {
        ok = FALSE;
        sprintf(buf, "element too big (array=%d)", i);
        DDD_PrintError('E', 9902, RegisterError(desc, 0, buf));
      }
#endif
    }
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
    else
    {
      if (CPP_STRUCT(desc) CPP_AND (e1->offset+e1->size > desc->size))
      {
        ok = FALSE;
        sprintf(buf, "element too big (offset=%d)", e1->offset);
        DDD_PrintError('E', 9903, RegisterError(desc, 0, buf));
      }
    }
#endif
  }
  return ok;
}



/*
        constructor for ELEM_DESC
 */

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
static void ConstructEl (ELEM_DESC *elem, int t, int o, size_t s, DDD_TYPE rt)
{
  elem->type    = t;
  elem->offset  = o;
  elem->size    = s;

  /*
          for OBJPTR elements, store referenced DDD_TYPE here.
          the default is EL_DDDHDR, i.e., if this feature is
          not used, the DDD_HDR will be assumed to be at the
          beginning of each structure (offsetHeader==0).
   */
  EDESC_SET_REFTYPE(elem, rt);
  elem->reftypeHandler = NULL;

  /*
          for GBITS elements, store array of bits. 1=GDATA,
          0=LDATA.
   */
  if (t==EL_GBITS)
  {
    elem->gbits = (unsigned char *) AllocFix(s);
    if (elem->gbits==NULL)
    {
      DDD_PrintError('E', 9932, STR_NOMEM " for EL_GBITS array");
      HARD_EXIT;
    }
  }
}
#endif


#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
static void ConstructEl (ELEM_DESC *elem, int t, char *a, size_t s, DDD_TYPE rt)
{
  elem->type    = t;
  elem->size    = s;

#ifdef CPP_FRONTEND
  /* in CPP, the first array entry is referenced by index 0. */
  elem->array   = a;
#else
  /* the size of one array entry (=s) is subtracted here,
     because we want to reference the first array entry
     by index 1. this is due to F77 conventions, where
     the first entry always is 1. in C, we would use 0
     as first index.
   */
  elem->array   = a - s;
#endif


  /*
          for OBJPTR elements, store referenced DDD_TYPE here.
          the default is EL_DDDHDR, i.e., if this feature is
          not used, the DDD_HDR will be assumed to be at the
          beginning of each structure (offsetHeader==0).
   */
  EDESC_SET_REFTYPE(elem, rt);
  elem->reftypeHandler = NULL;

  if (t==EL_GBITS)
  {
    /* TODO: GBITS could be supported also for F_FRONTEND and
            CPP_FRONTEND with STORAGE_ARRAY */
    DDD_PrintError('E', 9931, "EL_GBITS currently not supported");
    HARD_EXIT;
  }
}
#endif




/*
        register previously defined TYPE_DESC during TypeDefine
 */

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
static int RecursiveRegister (TYPE_DESC *desc,
                              int i, DDD_TYPE typ, int offs, int argno)
{
  TYPE_DESC *d2 = &(theTypeDefs[typ]);
  int j;
  char       *errtxt;

  /* inherit elements of other ddd-type */
  for(j=0; j<d2->nElements && i<MAX_ELEMDESC; j++, i++)
  {
    ConstructEl(&desc->element[i],
                d2->element[j].type,
                d2->element[j].offset + offs,
                d2->element[j].size,
                EDESC_REFTYPE(&(d2->element[j])));
    if (CheckBounds(desc, &desc->element[i], argno) == ERROR)
      return(ERROR);
  }

  /* inherit other properties */
  desc->nPointers += d2->nPointers;
  if (d2->hasHeader)
  {
    if (!desc->hasHeader)
    {
      desc->hasHeader = TRUE;
      desc->offsetHeader = d2->offsetHeader + offs;
    }
    else
    {
      if (desc->offsetHeader == d2->offsetHeader+offs)
      {
        errtxt=RegisterError(desc,argno, "two DDD_HDRs, same offset");
        DDD_PrintError('W', 9904, errtxt);
      }
      else
      {
        errtxt=RegisterError(desc,argno, "only one DDD_HDR allowed");
        DDD_PrintError('E', 9905, errtxt);
        return(ERROR);
      }
    }
  }

  return i;
}
#endif


#if defined(CPP_FRONTEND)
static int RecursiveRegister (TYPE_DESC *desc,
                              int i, DDD_TYPE typ, char *adr, int argno)
{
  TYPE_DESC *d2 = &(theTypeDefs[typ]);
  char       *errtxt;

  if (CPP_ARRAY(d2))
  {
    errtxt=RegisterError(desc,argno,
                         "cannot include array-like type into array-like type");
    DDD_PrintError('W', 9952, errtxt);
    return(ERROR);
  }

  ConstructEl(&desc->element[i], typ, adr, d2->size, 0);
  if (CheckBounds(desc, &desc->element[i], argno) == ERROR)
    return(ERROR);

  desc->size += d2->size;

  /* inherit other properties */
  desc->nPointers += d2->nPointers;
  if (d2->hasHeader)
  {
    if (!desc->hasHeader)
    {
      desc->hasHeader = TRUE;
      desc->elemHeader = i;
      desc->offsetHeader = d2->offsetHeader;
    }
    else
    {
      errtxt=RegisterError(desc,argno, "only one DDD_HDR allowed");
      DDD_PrintError('E', 9953, errtxt);
      return(ERROR);
    }
  }

  return i+1;
}
#endif



/*
        constructor for TYPE_DESC
 */

static void ConstructDesc (TYPE_DESC *desc)
{
  InitHandlers(desc);

  desc->nPointers = 0;
  desc->nElements = 0;
  desc->cmask     = NULL;
#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
  desc->size = 0;
#endif
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  desc->hasHeader = FALSE;
  desc->offsetHeader = 0;
#else
  if (! (desc->hdr = AllocHdr(sizeof(DDD_HEADER) * desc->arraySize) ) )
  {
    DDD_PrintError('E', 9999,
                   RegisterError(desc,0, STR_NOMEM));
    HARD_EXIT;             /*return;*/
  }
#endif
}



#ifndef DebugNoStructCompress

/*
        remove first element in given ELEM_DESC-list,
        adjust remaining elements
 */

static void DeleteFirstEl (ELEM_DESC *elarray, int n)
{
  int i;

  for(i=1; i<n; i++)
  {
    elarray[i-1] = elarray[i];
  }
}

#endif



/*
        normalize ELEM_DESC-list of given TYPE_DESC

        this consists of two parts:
                1) sort ELEM_DESC-list by offset (necessary!!)
                2) compress ELEM_DESC-list according to a set of rules
 */

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)

static int NormalizeDesc (TYPE_DESC *desc)
{
  ELEM_DESC  *elems = desc->element;

  /* sort element array by offset */
  qsort(elems, desc->nElements, sizeof(ELEM_DESC), sort_el_offset);

  /* check for overlapping elements */
  if (! CheckOverlapEls(desc))
    return FALSE;


#       ifndef DebugNoStructCompress
  {
    /* compile this only if Debug-flag not set */
    int i;

    /* compress element description */
    for(i=0; i<desc->nElements-1; i++)
    {
      size_t realsize;

      /* 1) type must be equal */
      if (elems[i].type != elems[i+1].type)
        continue;

      /* 2) nothing can melt into DDD_HEADER */
      if (desc->hasHeader && elems[i+1].offset==desc->offsetHeader)
        continue;

      /* 3) gap between elements is allowed only for EL_LDATA */
      if ((elems[i].offset+elems[i].size != elems[i+1].offset) &&
          (elems[i].type!=EL_LDATA) )
        continue;

      /* 4) EL_OBJPTRs with different reftypes can't be compressed */
      if (elems[i].type==EL_OBJPTR &&
          ((EDESC_REFTYPE(elems+i)!= EDESC_REFTYPE(elems+(i+1))) ||
           (EDESC_REFTYPE(elems+i)==DDD_TYPE_BY_HANDLER)) )
        continue;

      /* 5) EL_GBITS cant be compressed */
      if (elems[i].type == EL_GBITS)
        continue;

      /* all conditions fit: compress elements */
      realsize = elems[i+1].offset - elems[i].offset;
      elems[i].size = realsize + elems[i+1].size;

      desc->nElements--;
      DeleteFirstEl(&elems[i+1], desc->nElements - i);

      i--;                   /* skip one element back and try again */
    }
  }
#       endif

  return TRUE;
}


/*
        compute copy-mask (for efficiency) and attach it to TYPE_DESC
 */

static void AttachMask (TYPE_DESC *desc)
{
  int i, k;
  ELEM_DESC *e;
  unsigned char  *mp;
  unsigned char mask;

  /* get storage for mask */
  desc->cmask = (unsigned char *)AllocFix(desc->size);
  if (desc->cmask==0)
  {
    DDD_PrintError('E', 9906,
                   RegisterError(desc,0, STR_NOMEM));
    HARD_EXIT;             /*return;*/
  }

  /* set default: EL_LDATA for unspecified regions (gaps) */
  for(i=0; i<desc->size; i++)
  {
    desc->cmask[i] = 0x00;                    /* dont-copy-flag */
  }

  /* create mask from element list */
  for(i=0; i<desc->nElements; i++)
  {
    e = &desc->element[i];
    mp = desc->cmask + e->offset;

    switch (e->type)
    {
    case EL_LDATA :
    case EL_OBJPTR :                        /* OBJPTR are LDATA!!! */
      mask = 0x00;                              /* dont-copy-flag */
      break;

    case EL_GDATA :
    case EL_DATAPTR :
      mask = 0xff;                              /* copy-flag */
      break;
    }

    for(k=0; k<e->size; k++)
    {
      if (e->type==EL_GBITS)
      {
        mp[k] = e->gbits[k];                           /* copy bitwise */
      }
      else
      {
        mp[k] = mask;
      }
    }
  }

#       ifdef DebugCopyMask
  if (me==master)
  {
    char buf[8];

    sprintf(cBuffer, "%4d: AttachMask for %s:", me, desc->name);
    for(i=0; i<desc->size; i++)
    {
      if (i%8==0)
      {
        strcat(cBuffer,"\n");
        DDD_PrintLine(cBuffer);
        sprintf(cBuffer,"  %4d:  ", i);
      }
      sprintf(buf, "%02x ", desc->cmask[i]);
      strcat(cBuffer, buf);
    }
    strcat(cBuffer,"\n");
    DDD_PrintLine(cBuffer);
  }
#       endif
}

#endif


#ifdef F_FRONTEND


void ClearHeaders (TYPE_DESC *desc)

{
  DDD_HDR hdr;
  int i;

  for (i = 0, hdr = desc->hdr; i < desc->arraySize; i++, hdr++)
    MarkHdrInvalid (hdr);

  /* skip entry 0, because F77-arrays start with index 1.
     index 0 will be used as null-pointer. */
  desc->nextFree = 1;
}

#endif



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_TypeDefine                                                */
/*                                                                          */
/* Purpose:   define object structure at runtime                            */
/*                                                                          */
/* Input:     typ, adr, t0, p0, s0, [r0 [, rh0]], t1, p1, s1 ...                    */
/*            with typ:  previously declared DDD_TYPE                       */
/*                 adr:  example struct address                             */
/*                 t:    element type                                       */
/*                 p:    pointer into example structure                     */
/*                 s:    size of element in byte                            */
/*                 r:    (only for EL_OBJPTR): referenced DDD_TYPE          */
/*                 rt:   (only for EL_OBJPTR and DDD_TYPE_BY_HANDLER):      */
/*                          handler for getting reftype on-the-fly.         */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
void DDD_TypeDefine (DDD_TYPE typ, ...)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::TypeDefine (DDD_TYPE typ, ...)
#endif
#ifdef F_FRONTEND
void DDD_TypeDefine (DDD_TYPE *ftyp, ...)
#endif

{
  TYPE_DESC *desc;
  size_t argsize;
  char      *argp;
  int argtyp, argno;
  DDD_TYPE argrefs;
  int i, nPtr;
  char      *errtxt;
  va_list ap;
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  char      *adr;
  char      *gbits;
#endif
#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
  int size;
  int offset;
#endif
#ifdef F_FRONTEND
  DDD_TYPE typ = *ftyp;
#endif

  /* TODO: only master should be able to define types, other
          procs should receive the correct definition from master.
          (with the current implementation inconsistencies might occur)
   */

  /* test whether typ is valid */
  if (typ>=nDescr)
  {
    DDD_PrintError('E', 9907,
                   "invalid DDD_TYPE in DDD_TypeDefine");
    HARD_EXIT;             /*return;*/
  }

  /* get object description */
  desc = &(theTypeDefs[typ]);
  desc->currTypeDefCall++;

  if (desc->mode!=DDD_TYPE_DECLARED && desc->mode!=DDD_TYPE_CONTDEF)
  {
    if (desc->mode==DDD_TYPE_DEFINED)
    {
      DDD_PrintError('E', 9908,
                     RegisterError(desc, 0, "DDD_TYPE already defined"));
    }
    else
    {
      DDD_PrintError('E', 9908,
                     RegisterError(desc, 0, "undeclared DDD_TYPE"));
    }
    HARD_EXIT;             /*return;*/
  }


  /* initialize TYPE_DESC struct, only on first call */
  if (desc->currTypeDefCall==1)
    ConstructDesc(desc);


#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  if (typ==0)        /* i.e. typ==EL_DDDHDR */
  {
    /* DDD_HDR also contains a DDD_HDR (sic!) */
    desc->hasHeader = TRUE;
  }

#       ifdef DebugTypeDefine
  sprintf(cBuffer,"   DDD_TypeDefine(%s/%d)\n",
          desc->name, desc->currTypeDefCall);
  DDD_PrintDebug(cBuffer);
#       endif


  /* start variable arguments after "typ"-parameter */
  va_start(ap, typ);
#endif

#ifdef C_FRONTEND
  adr = va_arg(ap, char *);
  argno = 2;
#endif

#ifdef CPP_FRONTEND
  if (CPP_STRUCT(desc))
  {
    adr = va_arg(ap, char *);
    argno = 2;
  }
  else
  {
    // STORAGE_ARRAY mode
    size   = 0;
    offset = 0;
    argno  = 1;
  }
#endif

#ifdef F_FRONTEND
  va_start(ap, ftyp);

  size   = 0;
  offset = sizeof(DDD_HEADER);
  argno  = 1;
#endif


  /* loop over variable argument list */

  i = desc->nElements;
  while ((i<MAX_ELEMDESC) &&
         ((argtyp= FTYPE va_arg(ap, int FTYPE))!=EL_END) && (argtyp!=EL_CONTINUE))
  {
    HandlerGetRefType arg_rt_handler;

    /* get the pointer to the object (no special treatment for fortran	*/
    /* needed )															*/
    argp = va_arg(ap, char *);
    argno+=2;

    /* handle several types of ELEM_DESCs */
    switch (argtyp)
    {
    /* 1) ELEM_DESC is a pointer [array] */
    case EL_OBJPTR :
    case EL_DATAPTR :
      /* get third argument of this definition line */
      argsize = FTYPE va_arg(ap, size_t FTYPE); argno++;

      /* EL_OBJPTR have to be specified with target type */
      if (argtyp==EL_OBJPTR)
      {

        /* get fourth argument: referenced DDD_TYPE */
        argrefs = FTYPE va_arg(ap, DDD_TYPE FTYPE); argno++;

                                        #ifdef C_FRONTEND
        /* check whether target type is DDD_TYPE_BY_HANDLER */
        /* this is currently supported only by C_FRONTEND */
        if (argrefs==DDD_TYPE_BY_HANDLER)
        {
          arg_rt_handler = va_arg(ap, HandlerGetRefType);
          argno++;
        }
        else
                                        #endif
        {
          /* check whether target type is valid */
          if (argrefs>=nDescr ||
              theTypeDefs[argrefs].mode==DDD_TYPE_INVALID)
          {
            errtxt=RegisterError(desc,argno,
                                 "referencing invalid DDD_TYPE");
            DDD_PrintError('E', 9909, errtxt);
            HARD_EXIT;                                             /*return;*/
          }
        }
      }
      else
      {
        /* to target type for EL_DATAPTR */
        argrefs = EL_DDDHDR;
      }

      /* compute #pointers (in array) */
      nPtr = argsize / sizeof(void *);

      /* check for plausibility */
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      if (nPtr*sizeof(void *) != argsize)
#endif
#ifdef F_FRONTEND
      if (sizeof (DDD_OBJ) != argsize)
#endif
      {
        errtxt=RegisterError(desc,argno, "invalid sizeof");
        DDD_PrintError('E', 9910, errtxt);
        HARD_EXIT;                                 /*return;*/
      }

      /* remember #pointers */
      desc->nPointers += nPtr;


      /* initialize ELEM_DESC */
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      ConstructEl(&desc->element[i],
                  argtyp, (int)(argp-adr), argsize, argrefs);

      /* set reftype-handler function pointer, if any */
      if (argrefs==DDD_TYPE_BY_HANDLER)
        desc->element[i].reftypeHandler = arg_rt_handler;
#endif
#ifdef F_FRONTEND
      size += argsize;
      ConstructEl(&desc->element[i],
                  argtyp, (int)argp, argsize, argrefs);

      desc->element[i].msgoffset = offset;
      offset += argsize;
#endif
      if (CheckBounds(desc, &desc->element[i], argno) == ERROR)
        return;
      i++;

#                               ifdef DebugTypeDefine
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      sprintf(cBuffer,"    PTR, %05d, %06d\n",
              argp-adr, argsize);
#endif
#ifdef F_FRONTEND
      sprintf(cBuffer,"    PTR, %05d, %06d\n",
              argp, argsize);
#endif
      DDD_PrintDebug(cBuffer);
#                               endif

      break;


    /* 2) ELEM_DESC is global or local data */
    case EL_GDATA :
    case EL_LDATA :
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      /* get third argument of this definition line */
      argsize = va_arg(ap, size_t); argno++;

      /* initialize ELEM_DESC */
      if (CPP_STRUCT(desc) CPP_AND TRUE)
      {
        ConstructEl(&desc->element[i],
                    argtyp, (int)(argp-adr), argsize, 0);
      }
#endif
#ifdef F_FRONTEND
      argsize = *(long *) va_arg(ap, size_t); argno++;
#endif
#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
      if (CPP_ARRAY(desc) CPP_AND TRUE)
      {
        size += argsize;
        ConstructEl(&desc->element[i],
                    argtyp, argp, argsize, 0);

        if (argtyp == EL_GDATA)
        {
          desc->element[i].msgoffset = offset;
          offset += argsize;
        }
      }
#endif
      if (CheckBounds(desc, &desc->element[i], argno) == ERROR)
        return;
      i++;

#                               ifdef DebugTypeDefine
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      if (CPP_STRUCT(desc) CPP_AND TRUE)
      {
        sprintf(cBuffer,"    DAT, %05d, %06d\n",
                argp-adr, argsize);
      }
#endif
#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
      if (CPP_ARRAY(desc) CPP_AND TRUE)
      {
        sprintf(cBuffer,"    DAT, %08x, %06d\n",
                argp, argsize);
      }
#endif
      DDD_PrintDebug(cBuffer);
#                               endif

      break;



    /* 3) ELEM_DESC is bitwise global or local */
    case EL_GBITS :
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      /* get third argument of this definition line */
      argsize = va_arg(ap, size_t); argno++;

      /* initialize ELEM_DESC */
      ConstructEl(&desc->element[i],
                  argtyp, (int)(argp-adr), argsize, 0);

      /* read forth arg from cmdline */
      gbits = va_arg(ap, char *); argno++;

      /* fill gbits array, read forth arg from cmdline */
      memcpy(desc->element[i].gbits, gbits, argsize);
#endif
#ifdef F_FRONTEND
      /* TODO */
      DDD_PrintError('E', 9930, "EL_GBITS not supported in F_FRONTEND");
      HARD_EXIT;                           /*return;*/
#endif
      if (CheckBounds(desc, &desc->element[i], argno) == ERROR)
        return;

#                               ifdef DebugTypeDefine
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      sprintf(cBuffer,"   BITS, %05d, %06d, ",
              argp-adr, argsize);
      {
        int ii;
        char buf[5];
        for(ii=0; ii<argsize; ii++)
        {
          sprintf(buf, "%02x ", (int)desc->element[i].gbits[ii]);
          strcat(cBuffer, buf);
        }
        strcat(cBuffer, "\n");
      }
      DDD_PrintDebug(cBuffer);
#endif
#                               endif

      i++;

      break;


    /* 4) ELEM_DESC is a recursively defined DDD_TYPE */
    default :
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      /* hierarchical registering of known ddd_types */
      /* no third argument here */

      /* check for plausibility of given DDD_TYPE */
      if (argtyp<0 || argtyp>=nDescr || argtyp==typ)
      {
        char buf[40];
        sprintf(buf,"undefined DDD_TYPE=%d", argtyp);
        errtxt=RegisterError(desc,argno-1,buf);
        DDD_PrintError('E', 9911, errtxt);
        HARD_EXIT;                                 /*return;*/
      }

      /* check whether given DDD_TYPE has been defined already */
      if (theTypeDefs[argtyp].mode==DDD_TYPE_DEFINED)
      {
        if (CPP_STRUCT(desc) CPP_AND TRUE)
        {
          /* do recursive TypeDefine */
          i = RecursiveRegister(desc,
                                i, argtyp, (int)(argp-adr), argno);
          if (i==ERROR) HARD_EXIT;                                       /* return; */

                                                #ifdef DebugTypeDefine
          sprintf(cBuffer,"    %3d, %05d, %06d\n",
                  argtyp, argp-adr, theTypeDefs[argtyp].size);
          DDD_PrintDebug(cBuffer);
                                                #endif
        }
#ifdef CPP_FRONTEND
        else
        {
          /* do recursive TypeDefine */
          i = RecursiveRegister(desc, i, argtyp, argp, argno);
          if (i==ERROR) HARD_EXIT;                                       /* return; */

                                                #ifdef DebugTypeDefine
          sprintf(cBuffer,"    %3d, %08x, %06d\n",
                  argtyp, argp, theTypeDefs[argtyp].size);
          DDD_PrintDebug(cBuffer);
                                                #endif
        }
#endif
      }
      else
      {
        char buf[40];
        sprintf(buf,"undefined DDD_TYPE %s",
                theTypeDefs[argtyp].name);
        errtxt=RegisterError(desc,argno-1,buf);
        DDD_PrintError('E', 9912, errtxt);
        HARD_EXIT;                                 /*return;*/
      }

#endif
#ifdef F_FRONTEND
      errtxt=RegisterError(desc,argno,"recursive DDD_TYPE not impl");
      DDD_PrintError('E', 9912, errtxt);
      HARD_EXIT;                           /*return;*/
#endif
      break;
    }
  }


  /* check whether loop has come to a correct end */
  if (i>=MAX_ELEMDESC && argtyp!=EL_END && argtyp!=EL_CONTINUE)
  {
    errtxt=RegisterError(desc,0, "too many elements");
    DDD_PrintError('E', 1150, errtxt);
    HARD_EXIT;             /*return;*/
  }


  /* remember #elements in TYPE_DESC */
  desc->nElements = i;


  if (argtyp==EL_END)        /* and not EL_CONTINUE */
  {
    /* compute aligned object length */
                #if defined(C_FRONTEND)
    desc->size = (size_t) (va_arg(ap, char *) - adr);
    desc->size = CEIL(desc->size);
                #endif
                #if defined(CPP_FRONTEND)
    if (CPP_STRUCT(desc))
    {
      desc->size = va_arg(ap, char *) - adr;
      desc->size = CEIL(desc->size);
    }
    else
    {
      desc->size += size;
    }
                #endif
                #if defined(F_FRONTEND)
    desc->size = size + (desc->hdr ? sizeof (DDD_HEADER) : 0);
                #endif


#       if defined(C_FRONTEND) || defined(CPP_FRONTEND)
    if (CPP_STRUCT(desc) CPP_AND TRUE)
    {
      /* do normalization */
      if (! NormalizeDesc(desc))
        HARD_EXIT;                                 /*return;*/

      /* attach copy-mask for efficient copying */
      AttachMask(desc);
    }
#               endif
#               ifdef F_FRONTEND
    ClearHeaders (desc);
#       endif

    /* change TYPE_DESC state to DEFINED */
    desc->mode = DDD_TYPE_DEFINED;
  }
  else        /* argtyp==EL_CONTINUE */
  {
    /* change TYPE_DESC state to CONTDEF */
    desc->mode = DDD_TYPE_CONTDEF;
  }

  va_end(ap);
}




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_TypeDeclare                                               */
/*                                                                          */
/* Purpose:   declare DDD_TYPE at runtime                                   */
/*                                                                          */
/* Input:     name:  object name string                                     */
/*                                                                          */
/* Output:    id of object-description                                      */
/*            -1 on error                                                   */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
DDD_TYPE DDD_TypeDeclare (char *name)
#endif
#ifdef CPP_FRONTEND
DDD_TYPE DDD_Library::TypeDeclareStruct (char *name)
#endif
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
{
  TYPE_DESC *desc = &(theTypeDefs[nDescr]);

  /* check whether there is one more DDD_TYPE */
  if (nDescr==MAX_TYPEDESC)
  {
    DDD_PrintError('E', 9913, "no more DDD_TYPEs in DDD_TypeDeclare()");
    HARD_EXIT;             /*return(ERROR);*/
  }

  /* set status to DECLARED and remember textual type name */
  desc->mode = DDD_TYPE_DECLARED;
  desc->name = name;

  desc->prioMatrix  = NULL;
  desc->prioDefault = PRIOMERGE_DEFAULT;

#ifdef CPP_FRONTEND
  desc->storage  = STORAGE_STRUCT;
#endif

  /* increase #DDD_TYPEs, but return previously defined one */
  nDescr++; return(nDescr-1);
}
#endif



#ifdef CPP_FRONTEND
DDD_TYPE DDD_Library::TypeDeclareIndex (int size, char *name)
#endif
#ifdef F_FRONTEND
void DDD_TypeDeclare (char *name, int *size, DDD_TYPE *type)
#endif
#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
{
  TYPE_DESC *desc = &(theTypeDefs[nDescr]);

  /* check whether there is one more DDD_TYPE */
  if (nDescr==MAX_TYPEDESC)
  {
#ifdef CPP_FRONTEND
    DDD_PrintError('E', 9913, "no more DDD_TYPEs in DDD_TypeDeclare()");
    HARD_EXIT;             /*return(ERROR);*/
#else
    *type = -1;
    return;
#endif
  }

  /* set status to DECLARED and remember textual type name */
  desc->mode = DDD_TYPE_DECLARED;
  desc->name = name;

  desc->prioMatrix  = NULL;
  desc->prioDefault = PRIOMERGE_DEFAULT;


#ifdef CPP_FRONTEND
  desc->storage  = STORAGE_ARRAY;
  desc->arraySize = size;

  /* increase #DDD_TYPEs, but return previously defined one */
  nDescr++; return(nDescr-1);
#else
  desc->arraySize = *size;

  *type = nDescr++;
  return;
#endif
}
#endif


#ifdef CPP_FRONTEND
void DDD_Library::TypeChangeName (DDD_TYPE id, char *name)
{
  /* check for plausibility */
  if (id>=nDescr)
  {
    sprintf(cBuffer, "invalid DDD_TYPE %d in DDD_TypeChangeName", id);
    DDD_PrintError('E', 9933, cBuffer);
    HARD_EXIT;             /*return;*/
  }

  theTypeDefs[id].name = name;
}
#endif


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_TypeDisplay                                               */
/*                                                                          */
/* Purpose:   show defined DDD_TYPE                                         */
/*                                                                          */
/* Input:     id: DDD_TYPE which will be shown                              */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
void DDD_TypeDisplay (DDD_TYPE id)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::TypeDisplay (DDD_TYPE id)
#endif
#ifdef F_FRONTEND
void DDD_TypeDisplay (DDD_TYPE *idf)
#endif

{
  int i;
  TYPE_DESC *desc;
#ifdef F_FRONTEND
  DDD_TYPE id = *idf;
#endif

  /* only master should display DDD_TYPEs */
  if (me==master)
  {
    /* check for plausibility */
    if (id>=nDescr)
    {
      sprintf(cBuffer, "invalid DDD_TYPE %d in DDD_TypeDisplay", id);
      DDD_PrintError('E', 9914, cBuffer);
      HARD_EXIT;                   /*return;*/
    }

    desc = &(theTypeDefs[id]);
    if (desc->mode != DDD_TYPE_DEFINED)
    {
      sprintf(cBuffer, "undefined DDD_TYPE %d in DDD_TypeDisplay", id);
      DDD_PrintError('E', 9915, cBuffer);
      HARD_EXIT;                   /*return;*/
    }

    /* print header */
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
    if (CPP_STRUCT(desc) CPP_AND TRUE)
    {
      sprintf(cBuffer, "/ Structure of %s--object '%s', id %d, %d byte\n",
              desc->hasHeader ? "DDD" : "data",
              desc->name, id, desc->size);
    }
#endif
#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
    if (CPP_ARRAY(desc) CPP_AND TRUE)
    {
      sprintf(cBuffer,
              "/ Structure of %s--object '%s', id %d, %d byte, %d elements\n",
                                #ifdef F_FRONTEND
              desc->hdr ? "DDD" : "data",
                                #endif
                                #ifdef CPP_FRONTEND
              desc->hasHeader ? "DDD" : "data",
                                #endif
              desc->name, id, desc->size, desc->arraySize);
    }
#endif
    DDD_PrintLine(cBuffer);

    DDD_PrintLine(
      "|--------------------------------------------------------------\n");

    /* print one line for each element */
    for(i=0; i<desc->nElements; i++)
    {
      ELEM_DESC *e = &desc->element[i];
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      if (CPP_STRUCT(desc) CPP_AND TRUE)
      {
        int realnext = (i==desc->nElements-1) ? desc->size : (e+1)->offset;
        int estinext = e->offset+e->size;

        /* handle gap at the beginning */
        if (i==0 && e->offset!=0)
        {
          sprintf(cBuffer, "|%5d %5d    gap (local data)\n", 0, e->offset);
          DDD_PrintLine(cBuffer);
        }


        /* do visual compression of elems inherited from DDD_HDR */
        if (id==EL_DDDHDR ||
            (!desc->hasHeader) ||
            e->offset < desc->offsetHeader ||
            e->offset >= desc->offsetHeader+theTypeDefs[EL_DDDHDR].size)
        {
          sprintf(cBuffer, "|%5d %5d    ", e->offset, e->size);


          /* print one line according to type */
          switch (e->type)
          {
          case EL_GDATA : strcat(cBuffer, "global data\n"); break;
          case EL_LDATA : strcat(cBuffer, "local data\n"); break;
          case EL_DATAPTR : strcat(cBuffer, "data pointer\n"); break;
          case EL_OBJPTR :
            if (EDESC_REFTYPE(e)!=DDD_TYPE_BY_HANDLER)
            {
              sprintf(cBuffer, "%sobj pointer (refs %s)\n",
                      cBuffer,
                      theTypeDefs[EDESC_REFTYPE(e)].name);
            }
            else
            {
              sprintf(cBuffer,
                      "%sobj pointer (reftype on-the-fly)\n",
                      cBuffer);
            }
            break;
          case EL_GBITS : strcat(cBuffer, "bitwise global: ");
            {
              int ii;
              char buf[5];
              for(ii=0; ii<e->size; ii++)
              {
                sprintf(buf, "%02x ",
                        (int)e->gbits[ii]);
                strcat(cBuffer, buf);
              }
              strcat(cBuffer, "\n");
            }
            break;
          }
          DDD_PrintLine(cBuffer);


          /* handle gap */
          if (estinext != realnext)
          {
            sprintf(cBuffer, "|%5d %5d    gap (local data)\n",
                    estinext, realnext-estinext);
            DDD_PrintLine(cBuffer);
          }
        }
        else
        {
          /* handle included DDD_HDR */
          if (e->offset == desc->offsetHeader)
          {
            sprintf(cBuffer, "|%5d %5d    ddd-header\n",
                    e->offset, theTypeDefs[EL_DDDHDR].size);
            DDD_PrintLine(cBuffer);
          }
        }
      }
#endif
#if defined (F_FRONTEND) || defined(CPP_FRONTEND)
      if (CPP_ARRAY(desc) CPP_AND TRUE)
      {
        sprintf(cBuffer, "|%5d %5d    ", i, e->size);

        /* print one line according to type */
        switch (e->type)
        {
        case EL_GDATA : strcat(cBuffer, "global data\n"); break;
        case EL_LDATA : strcat(cBuffer, "local data\n"); break;
        case EL_DATAPTR : strcat(cBuffer, "data pointer\n"); break;
        case EL_OBJPTR :
          sprintf(cBuffer, "%sobj pointer (refs %s, offset %d)\n",
                  cBuffer,
                  theTypeDefs[e->reftype].name,e->msgoffset);
          break;
        default :
#ifdef F_FRONTEND
          sprintf(cBuffer, "%sunknown elemtype %d\n",
                  cBuffer, e->type);
#endif
#ifdef CPP_FRONTEND
          sprintf(cBuffer, "%srecursive type %s (typeId=%d)\n",
                  cBuffer,
                  theTypeDefs[e->type].name, e->type);
#endif
        }
        DDD_PrintLine(cBuffer);
      }
#endif
    }
    DDD_PrintLine(
      "\\--------------------------------------------------------------\n");
  }
}




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_HandlerRegister                                           */
/*                                                                          */
/* Purpose:   register handlers for given DDD_TYPE                          */
/*                                                                          */
/* Input:     typeId, handler_name, func_ptr, ... (last two repeated)      */
/*            with: typeId: DDD_TYPE for which handlers will be registered */
/*                  handler_name:  one of HANDLER_xxx                       */
/*                  func_ptr:      pointer to function implementing handler */
/*            the argument list must be finished with HANDLER_END.          */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

static void InitHandlers (TYPE_DESC *desc)
{
  /* set all handler functions to default (=none) */
  desc->handlerLDATACONSTRUCTOR = NULL;
  desc->handlerDESTRUCTOR = NULL;
  desc->handlerDELETE = NULL;
  desc->handlerUPDATE = NULL;
  desc->handlerOBJMKCONS = NULL;
  desc->handlerSETPRIORITY = NULL;
  desc->handlerXFERCOPY = NULL;
  desc->handlerXFERDELETE = NULL;
  desc->handlerXFERGATHER = NULL;
  desc->handlerXFERSCATTER = NULL;
  desc->handlerXFERGATHERX = NULL;
  desc->handlerXFERSCATTERX = NULL;
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  desc->handlerXFERCOPYMANIP = NULL;
#endif
#ifdef F_FRONTEND
  desc->handlerALLOCOBJ = NULL;
  desc->handlerFREEOBJ = NULL;
#endif
}


#define HDLR_NAME LDATACONSTRUCTOR
#include "handler.ct"

#define HDLR_NAME DESTRUCTOR
#include "handler.ct"

#define HDLR_NAME DELETE
#include "handler.ct"

#define HDLR_NAME UPDATE
#include "handler.ct"

#define HDLR_NAME OBJMKCONS
#include "handler.ct"

#define HDLR_NAME SETPRIORITY
#include "handler.ct"

#define HDLR_NAME XFERCOPY
#include "handler.ct"

#define HDLR_NAME XFERDELETE
#include "handler.ct"

#define HDLR_NAME XFERGATHER
#include "handler.ct"

#define HDLR_NAME XFERSCATTER
#include "handler.ct"

#define HDLR_NAME XFERGATHERX
#include "handler.ct"

#define HDLR_NAME XFERSCATTERX
#include "handler.ct"

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
#define HDLR_NAME XFERCOPYMANIP
#include "handler.ct"
#endif


#ifdef F_FRONTEND
#define HDLR_NAME ALLOCOBJ
#include "handler.ct"

#define HDLR_NAME FREEOBJ
#include "handler.ct"
#endif



/**
        Registration of handler functions.

        {\em This function is supported for downward compatibility only.
        (Use new DDD\_SetHandlerXXX-functions instead. Advantage: static type
        checking for handler functions).}

        This function registers a list of handlers for a certain DDD
        object type. The list may contain each handler type {\tt HANDLER\_}
        as listed in this manual. Usually, \funk{HandlerRegister} will be called
        once after initialization of the corresponding \cod{DDD\_TYPE}.

        It is possible to define certain sets of handler functions and switch them
        at runtime; \eg, a set for load migration and another set for distributed
        grid refinement can be defined. Before executing the actual task, the
        handler set must be registered.

        variable argumentlist...

   @param typeId  DDD type of object for which the handlers will be registered.
 */

#ifdef C_FRONTEND
void DDD_HandlerRegister (DDD_TYPE typeId, ...)
{
#endif
#ifdef F_FRONTEND
void DDD_HandlerRegister (DDD_TYPE *fid, ...)
{
  DDD_TYPE typeId = *fid;
#endif
#if defined(C_FRONTEND) || defined(F_FRONTEND)
TYPE_DESC *desc = &(theTypeDefs[typeId]);
int idx;
va_list ap;

OLDSTYLE("DDD_HandlerRegister() supported for downward compatibility only.");
OLDSTYLE("  (Use new DDD_SetHandlerXXX-functions instead.");
OLDSTYLE("   Advantage: static type checking for handler functions)");

if (desc->mode != DDD_TYPE_DEFINED)
{
  DDD_PrintError('E', 9916,
                 "undefined DDD_TYPE in DDD_HandlerRegister()");
  HARD_EXIT;               /*return;*/
}

/* read argument list, fill object structure definition */
        #ifdef C_FRONTEND
va_start(ap, typeId);
        #endif
        #ifdef F_FRONTEND
va_start(ap, fid);
        #endif

while ((idx = FTYPE va_arg(ap, int FTYPE)) != HANDLER_END)
{
  switch(idx)
  {
  case HANDLER_LDATACONSTRUCTOR :
    desc->handlerLDATACONSTRUCTOR =
      va_arg(ap, HandlerLDATACONSTRUCTOR);
    break;
  case HANDLER_DESTRUCTOR :
    desc->handlerDESTRUCTOR =
      va_arg(ap, HandlerDESTRUCTOR);
    break;
  case HANDLER_DELETE :
    desc->handlerDELETE =
      va_arg(ap, HandlerDELETE);
    break;
  case HANDLER_UPDATE :
    desc->handlerUPDATE =
      va_arg(ap, HandlerUPDATE);
    break;
  case HANDLER_OBJMKCONS :
    desc->handlerOBJMKCONS =
      va_arg(ap, HandlerOBJMKCONS);
    break;
  case HANDLER_SETPRIORITY :
    desc->handlerSETPRIORITY =
      va_arg(ap, HandlerSETPRIORITY);
    break;
  case HANDLER_XFERCOPY :
    desc->handlerXFERCOPY =
      va_arg(ap, HandlerXFERCOPY);
    break;
  case HANDLER_XFERDELETE :
    desc->handlerXFERDELETE =
      va_arg(ap, HandlerXFERDELETE);
    break;
  case HANDLER_XFERGATHER :
    desc->handlerXFERGATHER =
      va_arg(ap, HandlerXFERGATHER);
    break;
  case HANDLER_XFERSCATTER :
    desc->handlerXFERSCATTER =
      va_arg(ap, HandlerXFERSCATTER);
    break;
  case HANDLER_XFERGATHERX :
    desc->handlerXFERGATHERX =
      va_arg(ap, HandlerXFERGATHERX);
    break;
  case HANDLER_XFERSCATTERX :
    desc->handlerXFERSCATTERX =
      va_arg(ap, HandlerXFERSCATTERX);
    break;
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  case HANDLER_XFERCOPYMANIP :
    desc->handlerXFERCOPYMANIP =
      va_arg(ap, HandlerXFERCOPYMANIP);
    break;
#endif
#ifdef F_FRONTEND
  case HANDLER_ALLOCOBJ :
    desc->handlerALLOCOBJ =
      va_arg(ap, HandlerALLOCOBJ);
    break;
  case HANDLER_FREEOBJ :
    desc->handlerFREEOBJ =
      va_arg(ap, HandlerFREEOBJ);
    break;
#endif
  default :
    DDD_PrintError('E', 9917,
                   "undefined HandlerId in DDD_HandlerRegister()");
    HARD_EXIT;
  }
}

va_end(ap);
}

#endif


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_InfoTypes                                                 */
/*                                                                          */
/* Purpose:   returns number of declared types                              */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    number of declared types                                      */
/*                                                                          */
/****************************************************************************/

int DDD_InfoTypes (void)
{
  return nDescr;
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_InfoHdrOffset                                             */
/*                                                                          */
/* Purpose:   returns offset of DDD_HEADER for a given DDD_TYPE in bytes    */
/*            NOTE: output will be invalid for DDD_TYPEs without header!    */
/*                                                                          */
/* Input:     typeId:  DDD_TYPE for which header offset is returned        */
/*                                                                          */
/* Output:    nheader offset for given type                                 */
/*                                                                          */
/****************************************************************************/

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)

int DDD_InfoHdrOffset (DDD_TYPE typeId)
{
  TYPE_DESC *desc = &(theTypeDefs[typeId]);

  return desc->offsetHeader;
}

#endif

/****************************************************************************/
/*                                                                          */
/* Function:  ddd_TypeMgrInit                                               */
/*                                                                          */
/* Purpose:   init TypeMgr module                                           */
/*            (DDD_HEADER is declared and defined as first DDD_TYPE)        */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#if defined(C_FRONTEND) || defined(F_FRONTEND)
void ddd_TypeMgrInit (void)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::ddd_TypeMgrInit (void)
#endif
{
  int i;

  /* set all theTypeDefs to INVALID, i.e., no DDD_TYPE has been defined */
  for(i=0; i<MAX_TYPEDESC; i++)
  {
    theTypeDefs[i].mode = DDD_TYPE_INVALID;
    theTypeDefs[i].currTypeDefCall = 0;
  }


  /* reset declared types */
  nDescr = 0;


#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  /* init DDD_HEADER as first type, with DDD_TYPE=0 */
  {
    DDD_TYPE hdr_type;

                #ifdef C_FRONTEND
    DDD_HEADER *hdr = 0;

    /* hdr_type will be EL_DDDHDR (=0) per default */
    hdr_type = DDD_TypeDeclare("DDD_HDR");
    DDD_TypeDefine(hdr_type, hdr,
                #else
    DDD_Object *obj = 0;
    DDD_HEADER *hdr = &(obj->_hdr);

                        /* hdr_type will be EL_DDDHDR (=0) per default */
                        hdr_type = TypeDeclare("DDD_Object");
                        TypeDefine(hdr_type, obj,
                #endif

                   EL_GDATA, &hdr->typ,     sizeof(hdr->typ),
                   EL_LDATA, &hdr->prio,    sizeof(hdr->prio),
                   EL_GDATA, &hdr->attr,    sizeof(hdr->attr),
                   EL_GDATA, &hdr->flags,   sizeof(hdr->flags),
                   EL_LDATA, &hdr->myIndex, sizeof(hdr->myIndex),
                   EL_GDATA, &hdr->gid,     sizeof(hdr->gid),

                #ifdef C_FRONTEND
                   EL_END,   hdr+1  );
                #else
                   EL_END,   obj+1  );
                #endif
  }
#endif
}


/****************************************************************************/
/*                                                                          */
/* Function:  ddd_TypeMgrExit                                               */
/*                                                                          */
/* Purpose:   exit and clean up TypeMgr module                              */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void ddd_TypeMgrExit (void)
{
  int i;

  /* free memory */
  for(i=0; i<nDescr; i++)
  {
    if (theTypeDefs[i].cmask!=NULL)
      FreeFix(theTypeDefs[i].cmask);
  }
}




/****************************************************************************/
