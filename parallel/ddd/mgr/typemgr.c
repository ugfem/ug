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

#include "dddi.h"


/* #define DebugTypeDefine */
/* #define DebugCopyMask */
/* #define DebugNoStructCompress */



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


/*
        sort pointers to ELEM_DESC according to their offset
 */

#ifdef C_FRONTEND
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
#ifdef C_FRONTEND
  if (el->offset<0)
  {
    DDD_PrintError('E', 9900,
                   RegisterError(desc,argno, "negative offset"));
    return(ERROR);
  }
#else
  if (!el->array)
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
#ifdef C_FRONTEND
      if (e1->offset+e1->size > e2->offset)
      {
        ok = FALSE;
        sprintf(buf, "element too big (offset=%d)", e1->offset);
        DDD_PrintError('E', 9902, RegisterError(desc, 0, buf));
      }
#else
      if (e1->array+(e1->size*desc->arraySize) > e2->array)
      {
        ok = FALSE;
        sprintf(buf, "element too big (array=%d)", i);
        DDD_PrintError('E', 9902, RegisterError(desc, 0, buf));
      }
#endif
    }
#ifdef C_FRONTEND
    else
    {
      if (e1->offset+e1->size > desc->size)
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

#ifdef C_FRONTEND
static void ConstructEl (ELEM_DESC *elem, int t, int o, size_t s, DDD_TYPE rt)
#else
static void ConstructEl (ELEM_DESC *elem, int t, char *a, size_t s, DDD_TYPE rt)
#endif
{
  elem->type    = t;
#ifdef C_FRONTEND
  elem->offset  = o;
#else
  /* the size of one array entry (=s) is subtracted here,
     because we want to reference the first array entry
     by index 1. this is due to F77 conventions, where
     the first entry always is 1. in C, we would use 0
     as first index.
   */
  elem->array   = a - s;
#endif
  elem->size    = s;

  /*
          for OBJPTR elements, store referenced DDD_TYPE here.
          the default is EL_DDDHDR, i.e., if this feature is
          not used, the DDD_HDR will be assumed to be at the
          beginning of each structure (offsetHeader==0).
   */
  elem->reftype = rt;
}




/*
        register previously defined TYPE_DESC during TypeDefine
 */

#ifdef C_FRONTEND
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
                d2->element[j].reftype);
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



/*
        constructor for TYPE_DESC
 */

static void ConstructDesc (TYPE_DESC *desc)
{
  int i;

  /* set all handler functions to default (=none) */
  for(i=0; i<HANDLER_MAX; i++)
    desc->handler[i] = 0;

  desc->nPointers = 0;
  desc->nElements = 0;
#ifdef C_FRONTEND
  desc->hasHeader = FALSE;
  desc->offsetHeader = 0;
#else
  if (! (desc->hdr = AllocHdr(sizeof(DDD_HEADER) * desc->arraySize) ) )
  {
    DDD_PrintError('E', 9999,
                   RegisterError(desc,0, "out of memory"));
    exit(1);             /*return;*/
  }
#endif
}



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



/*
        normalize ELEM_DESC-list of given TYPE_DESC

        this consists of two parts:
                1) sort ELEM_DESC-list by offset (necessary!!)
                2) compress ELEM_DESC-list according to a set of rules
 */

#ifdef C_FRONTEND

static int NormalizeDesc (TYPE_DESC *desc)
{
  ELEM_DESC  *elems = desc->element;
  int i;

  /* sort element array by offset */
  qsort(elems, desc->nElements, sizeof(ELEM_DESC), sort_el_offset);

  /* check for overlapping elements */
  if (! CheckOverlapEls(desc))
    return FALSE;


#       ifdef DebugNoStructCompress
  return TRUE;
#       endif

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
        (elems[i].reftype != elems[i+1].reftype))
      continue;

    /* all conditions fit: compress elements */
    realsize = elems[i+1].offset - elems[i].offset;
    elems[i].size = realsize + elems[i+1].size;

    desc->nElements--;
    DeleteFirstEl(&elems[i+1], desc->nElements - i);

    i--;             /* skip one element back and try again */
  }

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
                   RegisterError(desc,0, "out of memory"));
    exit(1);             /*return;*/
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
    case EL_OBJPTR :
      mask = 0x00;                              /* dont-copy-flag */
      break;

    case EL_GDATA :
    case EL_DATAPTR :
      mask = 0xff;                              /* copy-flag */
      break;
    }

    for(k=0; k<e->size; k++)
    {
      mp[k] = mask;
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

#else

#define MarkHdrInvalid(hdr)    OBJ_INDEX(hdr)=MAX_OBJ

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
/* Input:     typ, adr, t0, p0, s0, [r0], t1, p1, s1 ...                    */
/*            with typ:  previously declared DDD_TYPE                       */
/*                 adr:  example struct address                             */
/*                 t:    element type                                       */
/*                 p:    pointer into example structure                     */
/*                 s:    size of element in byte                            */
/*                 r:    (only for EL_OBJPTR): referenced DDD_TYPE          */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
void DDD_TypeDefine (DDD_TYPE typ, ...)
#else
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
#ifdef C_FRONTEND
  char      *adr;
#else
  DDD_TYPE typ = *ftyp;
  int size;
  int offset;
#endif

  /* TODO: only master should be able to define types, other
          procs should receive the correct definition from master.
          (with the current implementation inconsistencies might occur)
   */

  /* test whether typ is valid */
  if (typ<0 || typ>=nDescr)
  {
    DDD_PrintError('E', 9907,
                   "invalid DDD_TYPE in DDD_TypeDefine");
    exit(1);             /*return;*/
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
    exit(1);             /*return;*/
  }


  /* initialize TYPE_DESC struct, only on first call */
  if (desc->currTypeDefCall==1)
    ConstructDesc(desc);


#ifdef C_FRONTEND
  if (typ==0)        /* i.e. typ==EL_DDDHDR */
  {
    /* DDD_HDR also contains a DDD_HDR (!) */
    desc->hasHeader = TRUE;
  }

#       ifdef DebugTypeDefine
  sprintf(cBuffer,"   DDD_TypeDefine(%s/%d)\n",
          desc->name, desc->currTypeDefCall);
  DDD_PrintDebug(cBuffer);
#       endif


  /* start variable arguments after "typ"-parameter */
  va_start(ap, typ);

  adr = va_arg(ap, char *);
  argno = 2;
#else
  va_start(ap, ftyp);

  size = 0;
  argno = 1;
  offset = sizeof(DDD_HEADER);
#endif

  /* loop over variable argument list */
  i = desc->nElements;
#ifdef C_FRONTEND
  while ((i<MAX_ELEMDESC) &&
         ((argtyp=va_arg(ap, int))!=EL_END) && (argtyp!=EL_CONTINUE))
#else
  while ((i<MAX_ELEMDESC) &&
         ((argtyp=*va_arg(ap, int *))!=EL_END) && (argtyp!=EL_CONTINUE))
#endif
  {
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
#ifdef C_FRONTEND
      argsize = va_arg(ap, size_t); argno++;
#else
      argsize = *va_arg(ap, size_t *); argno++;
#endif
      /* EL_OBJPTR have to be specified with target type */
      if (argtyp==EL_OBJPTR)
      {
        /* get fourth argument: referenced DDD_TYPE */
#ifdef C_FRONTEND
        argrefs = va_arg(ap, DDD_TYPE); argno++;
#else
        argrefs = *va_arg(ap, DDD_TYPE *); argno++;
#endif
        /* check whether target type is valid */
        if (argrefs<0 || argrefs>=nDescr ||
            theTypeDefs[argrefs].mode==DDD_TYPE_INVALID)
        {
          errtxt=RegisterError(desc,argno,
                               "referencing invalid DDD_TYPE");
          DDD_PrintError('E', 9909, errtxt);
          exit(1);                                       /*return;*/
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
#ifdef C_FRONTEND
      if (nPtr*sizeof(void *) != argsize)
#else
      if (sizeof (DDD_OBJ) != argsize)
#endif
      {
        errtxt=RegisterError(desc,argno, "invalid sizeof");
        DDD_PrintError('E', 9910, errtxt);
        exit(1);                                 /*return;*/
      }

      /* remember #pointers */
      desc->nPointers += nPtr;


      /* initialize ELEM_DESC */
#ifdef C_FRONTEND
      ConstructEl(&desc->element[i],
                  argtyp, argp-adr, argsize, argrefs);
#else
      size += argsize;
      ConstructEl(&desc->element[i],
                  argtyp, argp, argsize, argrefs);

      desc->element[i].msgoffset = offset;
      offset += argsize;
#endif
      if (CheckBounds(desc, &desc->element[i], argno) == ERROR)
        return;
      i++;

#                               ifdef DebugTypeDefine
#ifdef C_FRONTEND
      sprintf(cBuffer,"    PTR, %05d, %06d\n",
              argp-adr, argsize);
#else
      sprintf(cBuffer,"    PTR, %05d, %06d\n",
              argp, argsize);
#endif
      DDD_PrintDebug(cBuffer);
#                               endif

      break;


    /* 2) ELEM_DESC is global or local data */
    case EL_GDATA :
    case EL_LDATA :
#ifdef C_FRONTEND
      /* get third argument of this definition line */
      argsize = va_arg(ap, size_t); argno++;

      /* initialize ELEM_DESC */
      ConstructEl(&desc->element[i],
                  argtyp, argp-adr, argsize, 0);
#else
      argsize = *(long *) va_arg(ap, size_t); argno++;

      size += argsize;
      ConstructEl(&desc->element[i],
                  argtyp, argp, argsize, 0);

      if (argtyp == EL_GDATA)
      {
        desc->element[i].msgoffset = offset;
        offset += argsize;
      }
#endif
      if (CheckBounds(desc, &desc->element[i], argno) == ERROR)
        return;
      i++;

#                               ifdef DebugTypeDefine
#ifdef C_FRONTEND
      sprintf(cBuffer,"    DAT, %05d, %06d\n",
              argp-adr, argsize);
#else
      sprintf(cBuffer,"    DAT, %05d, %06d\n",
              argp, argsize);
#endif
      DDD_PrintDebug(cBuffer);
#                               endif

      break;


    /* 3) ELEM_DESC is a recursively defined DDD_TYPE */
    default :
#ifdef C_FRONTEND
      /* hierarchical registering of known ddd_types */
      /* no third argument here */

      /* check for plausibility of given DDD_TYPE */
      if (argtyp<0 || argtyp>=nDescr || argtyp==typ)
      {
        char buf[40];
        sprintf(buf,"undefined DDD_TYPE=%d", argtyp);
        errtxt=RegisterError(desc,argno-1,buf);
        DDD_PrintError('E', 9911, errtxt);
        exit(1);                                 /*return;*/
      }

      /* check whether given DDD_TYPE has been defined already */
      if (theTypeDefs[argtyp].mode==DDD_TYPE_DEFINED)
      {
        /* do recursive TypeDefine */
        i = RecursiveRegister(desc, i, argtyp, argp-adr, argno);
        if (i==ERROR) exit(1);                                 /* return; */
      }
      else
      {
        char buf[40];
        sprintf(buf,"undefined DDD_TYPE %s",
                theTypeDefs[argtyp].name);
        errtxt=RegisterError(desc,argno-1,buf);
        DDD_PrintError('E', 9912, errtxt);
        exit(1);                                 /*return;*/
      }

#                               ifdef DebugTypeDefine
      sprintf(cBuffer,"    %3d, %05d, %06d\n",
              argtyp, argp-adr, theTypeDefs[argtyp].size);
      DDD_PrintDebug(cBuffer);
#                               endif

#else
      errtxt=RegisterError(desc,argno,"recursive DDD_TYPE");
      DDD_PrintError('E', 9912, errtxt);
      exit(1);                           /*return;*/
#endif
      break;
    }
  }


  /* check whether loop has come to a correct end */
  if (i>=MAX_ELEMDESC && argtyp!=EL_END && argtyp!=EL_CONTINUE)
  {
    errtxt=RegisterError(desc,0, "too many elements");
    DDD_PrintError('E', 1150, errtxt);
    exit(1);             /*return;*/
  }


  /* remember #elements in TYPE_DESC */
  desc->nElements = i;


  if (argtyp==EL_END)        /* and not EL_CONTINUE */
  {
    /* compute aligned object length */
#               ifdef C_FRONTEND
    desc->size = va_arg(ap, char *) - adr;
    desc->size = CEIL(desc->size);
                #else
    desc->size = size + (desc->hdr ? sizeof (DDD_HEADER) : 0);
#               endif


#               ifdef C_FRONTEND
    /* do normalization */
    if (! NormalizeDesc(desc))
      exit(1);                           /*return;*/

    /* attach copy-mask for efficient copying */
    AttachMask(desc);
#               else
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
#else
void DDD_TypeDeclare (char *name, int *size, DDD_TYPE *type)
#endif
{
  TYPE_DESC *desc = &(theTypeDefs[nDescr]);

  /* check whether there is one more DDD_TYPE */
  if (nDescr==MAX_TYPEDESC)
  {
#ifdef C_FRONTEND
    DDD_PrintError('E', 9913, "no more DDD_TYPEs in DDD_TypeDeclare()");
    exit(1);             /*return(ERROR);*/
#else
    *type = -1;
    return;
#endif
  }

  /* set status to DECLARED and remember textual type name */
  desc->mode = DDD_TYPE_DECLARED;
  desc->name = name;

  desc->prioMatrix = NULL;


#ifdef C_FRONTEND
  /* increase #DDD_TYPEs, but return previously defined one */
  nDescr++; return(nDescr-1);
#else
  desc->arraySize = *size;

  *type = nDescr++;
  return;
#endif
}




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
#else
void DDD_TypeDisplay (DDD_TYPE *idf)
#endif

{
  int i;
  TYPE_DESC *desc;
  char buff[50];
#ifdef F_FRONTEND
  DDD_TYPE id = *idf;
#endif

  /* only master should display DDD_TYPEs */
  if (me==master)
  {
    /* check for plausibility */
    if (id<0 || id>=nDescr)
    {
      sprintf(cBuffer, "invalid DDD_TYPE %d in DDD_TypeDisplay", id);
      DDD_PrintError('E', 9914, cBuffer);
      exit(1);                   /*return;*/
    }

    desc = &(theTypeDefs[id]);
    if (desc->mode != DDD_TYPE_DEFINED)
    {
      sprintf(cBuffer, "undefined DDD_TYPE %d in DDD_TypeDisplay", id);
      DDD_PrintError('E', 9915, cBuffer);
      exit(1);                   /*return;*/
    }

    /* print header */
#ifdef C_FRONTEND
    sprintf(cBuffer, "/ Structure of %s--object '%s', id %d, %d byte\n",
            desc->hasHeader ? "DDD" : "data",
            desc->name, id, desc->size);
#else
    sprintf(cBuffer, "/ Structure of %s--object '%s', id %d, %d byte, %d elemnts\n",
            desc->hdr ? "DDD" : "data",
            desc->name, id, desc->size, desc->arraySize);
#endif
    DDD_PrintLine(cBuffer);

    DDD_PrintLine("|----------------------------------------------------\n");

    /* print one line for each element */
    for(i=0; i<desc->nElements; i++)
    {
      ELEM_DESC *e = &desc->element[i];
#ifdef C_FRONTEND
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
          sprintf(cBuffer, "%sobj pointer (refs %s)\n",
                  cBuffer,
                  theTypeDefs[e->reftype].name);
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
#else
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
      }
      DDD_PrintLine(cBuffer);
#endif
    }
    DDD_PrintLine("\\----------------------------------------------------\n");
  }
}




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_HandlerRegister                                           */
/*                                                                          */
/* Purpose:   register handlers for given DDD_TYPE                          */
/*                                                                          */
/* Input:     type_id, handler_name, func_ptr, ... (last two repeated)      */
/*            with: type_id: DDD_TYPE for which handlers will be registered */
/*                  handler_name:  one of HANDLER_xxx                       */
/*                  func_ptr:      pointer to function implementing handler */
/*            the argument list must be finished with HANDLER_END.          */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
void DDD_HandlerRegister (DDD_TYPE type_id, ...)
#else
void DDD_HandlerRegister (DDD_TYPE *fid, ...)
#endif
{
#ifdef F_FRONTEND
  DDD_TYPE type_id = *fid;
#endif
  TYPE_DESC *desc = &(theTypeDefs[type_id]);
  int idx;
  va_list ap;

  if (desc->mode != DDD_TYPE_DEFINED)
  {
    DDD_PrintError('E', 9916,
                   "undefined DDD_TYPE in DDD_HandlerRegister()");
    exit(1);             /*return;*/
  }

  /* read argument list, fill object structure definition */
#ifdef C_FRONTEND
  va_start(ap, type_id);

  while ((idx=va_arg(ap, int))!=HANDLER_END)
    desc->handler[idx] = va_arg(ap, HandlerPtr);
#else
  va_start(ap, fid);

  while ((idx = *va_arg(ap, int *)) != HANDLER_END)
    desc->handler [idx] = va_arg(ap, HandlerPtr);
#endif

  va_end(ap);
}



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
/* Input:     type_id:  DDD_TYPE for which header offset is returned        */
/*                                                                          */
/* Output:    nheader offset for given type                                 */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND

int DDD_InfoHdrOffset (DDD_TYPE type_id)
{
  TYPE_DESC *desc = &(theTypeDefs[type_id]);

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

void ddd_TypeMgrInit (void)
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


#ifdef C_FRONTEND
  /* init DDD_HEADER as first type, with DDD_TYPE=0 */
  {
    DDD_HEADER *hdr;
    DDD_TYPE hdr_type;

    /* hdr_type will be EL_DDDHDR (=0) per default */
    hdr_type = DDD_TypeDeclare("DDD_HDR");
    DDD_TypeDefine(hdr_type, hdr,
                   EL_GDATA, &hdr->typ,     sizeof(hdr->typ),
                   EL_LDATA, &hdr->prio,    sizeof(hdr->prio),
                   EL_GDATA, &hdr->attr,    sizeof(hdr->attr),
                   EL_GDATA, &hdr->flags,   sizeof(hdr->flags),
                   EL_LDATA, &hdr->myIndex, sizeof(hdr->myIndex),
                   EL_GDATA, &hdr->gid,     sizeof(hdr->gid),
                   EL_END,   hdr+1
                   );
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
{}




/****************************************************************************/
