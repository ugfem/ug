// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      pfile.c		                                                */
/*                                                                          */
/* Purpose:   a nice utility for writing a single file from parallel		*/
/*            processes														*/
/*                                                                          */
/* Author:	  Peter Bastian                                                                                         */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*																			*/
/* History:   28.01.97    begin												*/
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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "compiler.h"
#include "general.h"
#include "fileopen.h"

#include "pfile.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#undef LOCAL_DEBUG

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

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* definition of functions													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   pfile - write single (text) file from parallel processes

   DESCRIPTION:
   This module provides some useful functions to write a single
   file from multiple parallel processes. It is assumed that the
   file is subdivided into several segments. Each segment consists
   of any number of items (strings). The items may be written either
   in a non-deterministic order or the items may be tagged with
   a key and are then written in lowest key first order. The deterministic
   and nondeterministic write function may be mixed within a segment :-).

   The end of each segment therefore is a point of global synchronization
   of the processes. The amount of data written by each process is not
   known in advance.

   The end of file treated like the end of a segment with a terminate
   flag raised by each process.

   The items are strings and the file is opened in ASCII mode.

   In sequential mode the functions are mapped to standard C library
   functions.

   The number of open files is determined by available memory and the operating
   system.

   pfile uses the Concentrate, Broadcast calls of ppif.

   FILE STRUCTURE:

   .vb
   segment 1
      item1 ... item n1
   segment 2
      item1 ... item n2
   ...
   segment m
      item1 ... item nm
   .ve

   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   pfile_open - open a parallel file

   SYNOPSIS:
   PFILE *pfile_open (char *name);

   PARAMETERS:
   .  name - filename, relative paths are allowed

   DESCRIPTION:
   Allocates a new PFILE data structure on all processors and
   returns a pointer to it. This function has to be called by
   all processes, but the file is only opened on the master.
   The PFILE data structure contains a buffer. The size of
   the buffer can be adjusted in pfile.h .

   RETURN VALUE:
   .n NULL if any errors are encountered, this is a global state
   .n valid pointer if allocate successful on all processors

   SEE ALSO:
   'pfile_sync', 'pfile_puts', 'pfile_close'
   D*/
/****************************************************************************/

PFILE *pfile_open (char *name)
{
        #ifndef ModelP
  return( (PFILE *) fileopen(name,"w") );
        #else
  PFILE *pf;
  FILE *stream;
  INT i,error;

  /* allocate PFILE data structure */
  error = 0;
  pf = malloc(sizeof(PFILE));
  if (pf == NULL) error=1;
  error = UG_GlobalMaxINT(error);
  if (error) {
    if (pf != NULL) free(pf);
    return(NULL);
  }

  /* open file */
  error = 0;
  if (me == master)
  {
    stream = fileopen(name,"w");
    if (stream == NULL) error = 1;
  }
  else
  {
    stream = NULL;
  }
  error = UG_GlobalMaxINT(error);
  if (error) {
    free(pf); return(NULL);
  }

  /* init structure */
  pf->stream = stream;

  pf->local.nchars = 0;
  pf->local.first_key = 0;
  pf->local.last_key = 0;

  for (i=0; i<PFILE_MAX_TREE; i++)
    pf->valid_state[i] = 0;

  pf->robin = 0;

  /* and return */
  return(pf);

        #endif
}


/****************************************************************************/
/*D
   pfile_master_puts - write string to file immediately

   SYNOPSIS:
   INT pfile_master_puts (PFILE *pf, char *s);

   PARAMETERS:
   .  pf - pointer to parallel file
   .  s - item to be written

   DESCRIPTION:
   Enables the master process to write a string to the
   output file immediately. This function should be called
   immediately after pfile_open or pfile_sync. Otherwise
   the position of the string withinm the segment is not
   predictable.

   RETURN VALUE:
   .n 0 if ok.
   .n >0 if any errors are encountered, this is a LOCAL state

   SEE ALSO:
   'pfile_sync', 'pfile_open', 'pfile_close'
   D*/
/****************************************************************************/

INT pfile_master_puts (PFILE *pf, char *s)
{
        #ifndef ModelP
  fputs(s,(FILE *) pf);
        #else
  if (me == master)
    fputs(s,pf->stream);
        #endif

  return(0);
}

#ifdef ModelP
/****************************************************************************/
/*
   flush_buffer - flush buffer in global sequence

   SYNOPSIS:
   static void flush_buffer (PFILE *pf)

   PARAMETERS:
   .  pf - pointer to parallel file

   DESCRIPTION:

   RETURN VALUE:
   .n 1 if all processors had local.nchars==0
   .n 0 else

   SEE ALSO:
   'pfile_sync', 'pfile_open', 'pfile_close'
 */
/****************************************************************************/
static INT flush_buffer (PFILE *pf)
{
  INT i,j,min_key,min_index,finish=0;
  char *buffer;
  INT nchars;

  /* copy local state to array */
  pf->state[degree] = pf->local;

  /* update downtree states */
  for (i=0; i<degree; i++)
    if (!pf->valid_state[i])
    {
      GetConcentrate(i,pf->state+i,sizeof(PFILE_STATE));
      pf->valid_state[i] = 1;
                        #ifdef LOCAL_DEBUG
      UserWriteF("receiving state from %d: n=%d, fst=%d, lst=%d\n",
                 i,pf->state[i].nchars,pf->state[i].first_key,pf->state[i].last_key);
                        #endif
    }
  /* now we have a valid state from all downtree nodes	*/
  /* and ourself.											*/

  /* determine smallest key */
  min_key = PFILE_MAX_INT;
  for (j=0; j<=degree; j++)
  {
    i = (pf->robin+j) % (degree+1);             /* makes load balancing in case of equal keys */
    if ( (pf->state[i].nchars>0) && (pf->state[i].first_key<min_key) )
    {
      min_key = pf->state[i].first_key;
      min_index = i;
    }
  }
  pf->robin = (pf->robin+1) % (degree+1);
  /* now we know the smallest key and its owner */
        #ifdef LOCAL_DEBUG
  UserWriteF("min_index=%d, min_key=%d\n",min_index,min_key);
        #endif

  /* now the termination protocol:							*/
  /* if no smallest key has been found subtree including me	*/
  /* is finished, if me==master then all are finished			*/
  if (min_key == PFILE_MAX_INT)
  {
    if (me == master)
      finish = 1;                   /* now all are finished */
    else {
      Concentrate(&pf->local,sizeof(PFILE_STATE));                   /* has nchars==0 !  */
      GetSpread(&finish,sizeof(INT));                                        /* the global state */
    }

    /* now we have finish from uptree */
    if (finish)
    {
      /* all are finished, inform downtree nodes */
      for (i=0; i<degree; i++) {
        Spread(i,&finish,sizeof(INT));
        pf->valid_state[i] = 0;                               /* state is invalid */
      }
    }
                #ifdef LOCAL_DEBUG
    UserWriteF("finish=%d\n",finish);
                #endif

    /* now we are either finished or we will find min_key == PFILE_MAX_INT	*/
    /* when we enter next time                                                                                  */
    return(finish);
  }

  /* a smallest key has been found, so buffer is written to	*/
  /* output file or passed uptree.							*/
  /* first, lets get the data                                                           */
  if (min_index==degree)
  {
    /* it is our own data */
                #ifdef LOCAL_DEBUG
    UserWriteF("Writing my own buffer, nchars=%d\n",pf->local.nchars);
                #endif
    buffer = pf->buffer;
    nchars = pf->local.nchars;
    pf->local.nchars = 0;             /* buffer is empty again */
    pf->local.first_key = pf->local.last_key = 0;
  }
  else
  {
    /* fetch downtree data, we have only state ! */
    GetConcentrate(min_index,pf->buffer2,pf->state[min_index].nchars+1);
                #ifdef LOCAL_DEBUG
    UserWriteF("Writing buffer from %d, nchars=%d\n",min_index,pf->state[min_index].nchars);
                #endif
    buffer = pf->buffer2;
    nchars = pf->state[min_index].nchars;
    pf->valid_state[min_index] = 0;             /* get new state next time */
  }
  /* now process data */
  if (me==master) {
    if (nchars>0) fputs(buffer,pf->stream);
  }
  else {
    Concentrate(pf->state+min_index,sizeof(PFILE_STATE));             /* send state */
    Concentrate(buffer,nchars+1);                  /* send data  */
  }

  /* Note: While sending data uptree we don't expect to get	*/
  /* an answer from uptree!									*/
  return(0);       /* we are not finished yet */
}
#endif


/****************************************************************************/
/*D
   pfile_puts - write string to parallel file

   SYNOPSIS:
   INT pfile_puts (PFILE *pf, char *s);

   PARAMETERS:
   .  pf - pointer to parallel file
   .  s - item to be written

   DESCRIPTION:
   Writes a string to a previously opened parallel file.
   It is ensured that the string s is written as a whole to
   the output file. The string is appended to the output buffer
   without increasing the key.

   Any number of calls to pfile_puts may be issued by a single process.

   RETURN VALUE:
   .n 0 if ok.
   .n >0 if any errors are encountered, this is a LOCAL state

   SEE ALSO:
   'pfile_sync', 'pfile_open', 'pfile_close'
   D*/
/****************************************************************************/

INT pfile_puts (PFILE *pf, char *s)
{
        #ifndef ModelP
  fputs(s,(FILE *) pf);
        #else
  FILE *stream;
  INT error,n;

  /* do not change keys here */

  /* put into buffer if possible */
  n = strlen(s);
  if (n+pf->local.nchars<PFILE_BUFFER_SIZE)
  {
    strcpy(&(pf->buffer[pf->local.nchars]),s);
    pf->local.nchars += n;
  }
  else
  {
    if (n>=PFILE_BUFFER_SIZE)
    {
      PrintErrorMessage('E',"pfile_puts","string larger than buffer");
      return(1);
    }

    /* flush buffer */
    if (flush_buffer(pf))
    {
      PrintErrorMessage('E',"pfile_puts","unxepected finish");
      return(1);
    }

    /* and try again */
    strcpy(&(pf->buffer[pf->local.nchars]),s);
    pf->local.nchars += n;
  }
        #endif

  return(0);
}



/****************************************************************************/
/*D
   pfile_tagged_puts - write tagged string to parallel file

   SYNOPSIS:
   INT pfile_tagged_puts (PFILE *pf, char *s, INT key);

   PARAMETERS:
   .  pf - pointer to parallel file
   .  s - item to be written
   .  key - key for this item

   DESCRIPTION:
   Writes a tagged item to the output file. Tags should be globally
   unique and locally increasing within each segment. If tags
   are out of order, a warning is issued and the item is treated
   as in order.

   Any number of calls to pfile_tagged_puts may be issued by a single process.

   RETURN VALUE:
   .n 0 if ok.
   .n >0 if any errors are encountered, this is a LOCAL state

   SEE ALSO:
   'pfile_sync', 'pfile_open', 'pfile_close'
   D*/
/****************************************************************************/

#ifdef ModelP
static INT append_buffer (PFILE *pf, char *s, INT key)
{
  INT n;

        #ifdef LOCAL_DEBUG
  UserWriteF("appending %s with key %d, nchars=%d, lst=%d\n",s,key,
             pf->local.nchars,pf->local.last_key);
        #endif
  n = strlen(s);
  if (n+pf->local.nchars<PFILE_BUFFER_SIZE)
  {
    /* initialize keys if buffer is empty */
    if (pf->local.nchars == 0) pf->local.first_key = key;

    /* append to buffer */
    strcpy(&(pf->buffer[pf->local.nchars]),s);
    pf->local.nchars += n;

    /* update last key */
    pf->local.last_key = key;

    return(1);             /* append successful */
  }

  if (n>=PFILE_BUFFER_SIZE)
  {
    PrintErrorMessage('E',"pfile_puts","string larger than buffer");
    return(1);             /* treat as appended ! */
  }

  return(0);       /* could not append */
}
#endif

INT pfile_tagged_puts (PFILE *pf, char *s, INT key)
{
        #ifndef ModelP
  fputs(s,(FILE *) pf);
        #else
  /* check order of keys */
  if ( (pf->local.nchars > 0) && (key < pf->local.last_key) )
  {
    PrintErrorMessage('W',"pfile_tagged_puts","keys locally not ordered");
    key = pf->local.last_key+1;             /* treat as in order ! */
  }

  if ( (pf->local.nchars > 0) && (key > pf->local.last_key+1) )       /* not consecutive case */
    if (flush_buffer(pf))
    {
      PrintErrorMessage('E',"pfile_tagged_puts","unxepected finish");
      return(1);
    }

  if (!append_buffer(pf,s,key))
  {
    if (flush_buffer(pf))
    {
      PrintErrorMessage('E',"pfile_tagged_puts","unxepected finish");
      return(1);
    }
    append_buffer(pf,s,key);
  }
        #endif

  return(0);
}


/****************************************************************************/
/*D
   pfile_sync - indicate end of file segment

   SYNOPSIS:
   INT pfile_sync (PFILE *pf);

   PARAMETERS:
   .  pf - pointer to parallel file

   DESCRIPTION:
   At the end of a segment this function has to be called by
   all processes. The function returns when all process
   have reached the end the current segment. Then the segment
   counter is increased.

   RETURN VALUE:
   .n 0 if ok.
   .n >0 if any errors are encountered, this is a LOCAL state

   SEE ALSO:
   'pfile_puts', 'pfile_open', 'pfile_close'
   D*/
/****************************************************************************/

INT pfile_sync (PFILE *pf)
{
        #ifdef ModelP
  /* wait until all processors reach end of segment */
  while (!flush_buffer(pf)) ;
        #endif

  return(0);
}



/****************************************************************************/
/*D
   pfile_close - indicate end of file

   SYNOPSIS:
   INT pfile_close (PFILE *pf);

   PARAMETERS:
   .  pf - pointer to parallel file

   DESCRIPTION:
   At the end of the file this function has to be called by
   all processes. The function returns when all process
   have reached the end of file. Then the file is closed
   and buffer space is released.

   RETURN VALUE:
   .n 0 if ok.
   .n >0 if any errors are encountered, this is a LOCAL state

   SEE ALSO:
   'pfile_puts', 'pfile_open', 'pfile_sync'
   D*/
/****************************************************************************/

INT pfile_close (PFILE *pf)
{
        #ifdef ModelP
  /* wait until all processors reach end of segment */
  while (!flush_buffer(pf)) ;

  /* close output file */
  if (me == master) fclose(pf->stream);

  /* free data structure */
  free(pf);
        #else
  fclose( (FILE *) pf );
        #endif

  return(0);
}
