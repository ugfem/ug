// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      memutil.h                                                     */
/*                                                                          */
/* Purpose:   data types for buffer management                              */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   970422 kb  begin                                              */
/*                                                                          */
/* Remarks:   the data types in this file are implemented as                */
/*            'pseudo classes', as C++ isn't used in this implementation.   */
/*            the member function names start with the class names; and the */
/*            first parameter is the object instance.                       */
/*                                                                          */
/****************************************************************************/

/* RCS_ID
   $Header$
 */


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __DDD_MEMUTIL_H__
#define __DDD_MEMUTIL_H__



/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

/*
        class Buffer

        a Buffer is a contiguous chunk of memory. it has two access levels:
        lower level:  alloc/free, same as malloc()/free(), mapped to
                      AllocMsg/FreeMsg.
        higher level: allocated size != used size, the buffer size can be decreased
                      (fast, but wastes memory) or increased (with alloc) after
                      allocation of the buffer.
 */


typedef struct _Buffer
{
  char    *buf;              /* pointer to memory for buffer                    */
  size_t size;               /* size of memory chunk                            */
  size_t used;               /* size of used memory inside memory, used<=size!  */
} Buffer;



/* macros (member functions) */

/* initialize empty buffer */
#define BufferInit(b)     { (b).buf=NULL; (b).size=(b).used=0; }

/* release buffer's memory, set buffer to empty */
#define BufferFree(b)     { if ((b).buf!=NULL)   \
                            { FreeMsg((b).buf,(b).size); BufferInit(b); }  }

/* allocate memory, doesn't control if there already is allocated memory */
/* note: s==0 isn't checked, AllocMsg==0 isn't check. TODO. */
#define BufferAlloc(b,s)  { (b).buf=(char *)AllocMsg(s); (b).size=(b).used=(s); }

/* empty buffer virtually */
#define BufferReset(b)    { (b).used=0; }

/* reuse buffer or increase size, if existing buffer space is too small */
#define BufferCreate(b,s) {                                   \
    if ((s)<=(b).size)                                \
      (b).used = (s);                               \
    else {                                            \
      if ((b).buf!=NULL) FreeMsg((b).buf,(b).size); \
      BufferAlloc((b),(s));                         \
    }  }

/* get pointer to buffer memory */
#define BufferMem(b)      ((b).buf)

/* get length of used part of buffer */
#define BufferLen(b)      ((b).used)

/* return true if buffer is empty */
#define BufferIsEmpty(b)  ((b).used==0)



/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/


/* memutil.c */


#endif
