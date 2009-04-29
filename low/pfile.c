// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      pfile.c                                                       */
/*                                                                          */
/* Purpose:   a nice utility for writing a single file from parallel        */
/*            processes                                                     */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*                                                                          */
/* History:   28.01.97    begin                                             */
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

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "compiler.h"
#include "ugdevices.h"
#include "general.h"
#include "fileopen.h"

#ifdef ModelP
#include "pargm.h"
USING_PPIF_NAMESPACE
#endif

#include "pfile.h"

USING_UG_NAMESPACE

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
/*                                                                          */
/* definition of functions													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/** \file
    \brief write single (text) file from parallel processes

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

   \section File Structure

   \verbatim
   segment 1
      item1 ... item n1
   segment 2
      item1 ... item n2
   ...
   segment m
      item1 ... item nm
   \endverbatim

 */
/****************************************************************************/

/****************************************************************************/
/** \brief Open a parallel file

   \param name - filename, relative paths are allowed

   Allocates a new PFILE data structure on all processors and
   returns a pointer to it. This function has to be called by
   all processes, but the file is only opened on the master.
   The PFILE data structure contains a buffer. The size of
   the buffer can be adjusted in pfile.h .

   \return <ul>
   <li> NULL if any errors are encountered, this is a global state </li>
   <li> valid pointer if allocate successful on all processors </li>
   </ul>

   \sa
   'pfile_sync', 'pfile_puts', 'pfile_close'
 */
/****************************************************************************/

NS_DIM_PREFIX PFILE * NS_DIM_PREFIX pfile_open (char *name)
{
#ifndef ModelP
  return( (PFILE *) fileopen(name,"w") );
#else
  PFILE *pf;
  FILE *stream;
  INT i,error;

  /* allocate PFILE data structure */
  error = 0;
  pf = (PFILE*)malloc(sizeof(PFILE));
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
/** \brief Write string to file immediately

   \param pf - pointer to parallel file
   \param s - item to be written

   Enables the master process to write a string to the
   output file immediately. This function should be called
   immediately after pfile_open or pfile_sync. Otherwise
   the position of the string withinm the segment is not
   predictable.

   \return <ul>
   <li> 0 if ok </li>
   <li> >0 if any errors are encountered, this is a LOCAL state </li>
   </ul>

   \sa
   'pfile_sync', 'pfile_open', 'pfile_close'
 */
/****************************************************************************/

INT NS_DIM_PREFIX pfile_master_puts (PFILE *pf, char *s)
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
/** \brief Flush buffer in global sequence

   \param pf - pointer to parallel file

   \return <ul>
   <li> 1 if all processors had local.nchars==0 </li>
   <li> 0 else </li>
   </ul>

   \sa
   'pfile_sync', 'pfile_open', 'pfile_close'
 */
/****************************************************************************/
static INT flush_buffer (NS_DIM_PREFIX PFILE *pf)
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
      GetConcentrate(i,pf->state+i,sizeof(NS_DIM_PREFIX PFILE_STATE));
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
      Concentrate(&pf->local,sizeof(NS_DIM_PREFIX PFILE_STATE));                   /* has nchars==0 !  */
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
    Concentrate(pf->state+min_index,sizeof(NS_DIM_PREFIX PFILE_STATE));             /* send state */
    Concentrate(buffer,nchars+1);                  /* send data  */
  }

  /* Note: While sending data uptree we don't expect to get	*/
  /* an answer from uptree!									*/
  return(0);       /* we are not finished yet */
}
#endif


/****************************************************************************/
/** \brief Write string to parallel file

   \param pf - pointer to parallel file
   \param s - item to be written

   Writes a string to a previously opened parallel file.
   It is ensured that the string s is written as a whole to
   the output file. The string is appended to the output buffer
   without increasing the key.

   Any number of calls to pfile_puts may be issued by a single process.

   \return <ul>
   <li> 0 if ok </li>
   <li> >0 if any errors are encountered, this is a LOCAL state </li>
   </ul>

   \sa
   'pfile_sync', 'pfile_open', 'pfile_close'
 */
/****************************************************************************/

INT NS_DIM_PREFIX pfile_puts (PFILE *pf, char *s)
{
#ifndef ModelP
  fputs(s,(FILE *) pf);
#else
  INT n;

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



#ifdef ModelP
static INT append_buffer (NS_DIM_PREFIX PFILE *pf, char *s, INT key)
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

/****************************************************************************/
/** \brief Write tagged string to parallel file

   \param pf - pointer to parallel file
   \param s - item to be written
   \param key - key for this item

   Writes a tagged item to the output file. Tags should be globally
   unique and locally increasing within each segment. If tags
   are out of order, a warning is issued and the item is treated
   as in order.

   Any number of calls to pfile_tagged_puts may be issued by a single process.

   \return <ul>
   <li> 0 if ok </li>
   <li> >0 if any errors are encountered, this is a LOCAL state </li>
   </ul>

   \sa
   'pfile_sync', 'pfile_open', 'pfile_close'
 */
/****************************************************************************/

INT NS_DIM_PREFIX pfile_tagged_puts (PFILE *pf, char *s, INT key)
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
/** \brief Indicate end of file segment

   \param pf - pointer to parallel file

   At the end of a segment this function has to be called by
   all processes. The function returns when all process
   have reached the end the current segment. Then the segment
   counter is increased.

   \return <ul>
   <li> 0 if ok </li>
   <li> >0 if any errors are encountered, this is a LOCAL state </li>
   </ul>

   \sa
   'pfile_puts', 'pfile_open', 'pfile_close'
 */
/****************************************************************************/

INT NS_DIM_PREFIX pfile_sync (PFILE *pf)
{
#ifdef ModelP
  /* wait until all processors reach end of segment */
  while (!flush_buffer(pf)) ;
#endif

  return(0);
}



/****************************************************************************/
/** \brief Indicate end of file

   \param pf - pointer to parallel file

   At the end of the file this function has to be called by
   all processes. The function returns when all process
   have reached the end of file. Then the file is closed
   and buffer space is released.

   \return <ul>
   <li> 0 if ok </li>
   <li> >0 if any errors are encountered, this is a LOCAL state </li>
   </ul>

   \sa
   'pfile_puts', 'pfile_open', 'pfile_sync'
 */
/****************************************************************************/

INT NS_DIM_PREFIX pfile_close (PFILE *pf)
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

/****************************************************************************/
/** \brief Open a parallel binary file

   \param name - filename, relative paths are allowed

   Allocates a new PFILE_BIN data structure on all processors and
   returns a pointer to it. This function has to be called by
   all processes, but the file is only opened on the master.
   The PFILE_BIN data structure contains 2 buffers, one for integers
   and one for floats. The size of the buffer can be adjusted
   in pfile.h .

   \return <ul>
   <li> NULL if any errors are encountered, this is a global state </li>
   <li> valid pointer if allocate successful on all processors </li>
   </ul>

   \sa
   'pfile_sync_bin', 'pfile_close_bin'
 */
/****************************************************************************/

NS_DIM_PREFIX PFILE_BIN* NS_DIM_PREFIX pfile_open_bin (char *name)
{
#ifndef ModelP
  return( (PFILE_BIN *) fileopen(name,"wb") );
#else
  PFILE_BIN *pf;
  FILE *stream;
  INT i,error;

  /* allocate PFILE_BIN data structure */
  error = 0;
  pf = (PFILE_BIN*)malloc(sizeof(PFILE_BIN));
  if (pf == NULL) error=1;
  error = UG_GlobalMaxINT(error);
  if (error) {
    if (pf != NULL) free(pf);
    return(NULL);
  }

  /* open file */
  error = 0;
  if (me == master) {
    stream = fileopen(name,"w");
    if (stream == NULL) error = 1;
  } else
    stream = NULL;
  error = UG_GlobalMaxINT(error);
  if (error) {
    free(pf);
    return(NULL);
  }

  /* init structure */
  pf->stream = stream;

  pf->local.nints = 0;
  pf->local.nfloats = 0;
  pf->local.first_key = 0;
  pf->local.last_key = 0;

  for (i=0; i<PFILE_MAX_TREE; i++)
    pf->valid_state[i] = 0;

  pf->robin = 0;

  /* and return */
  return(pf);
#endif
}

#ifdef ModelP
/****************************************************************************/
/* change for BYTES not finished ! */
static INT flush_buffer_bin (NS_DIM_PREFIX PFILE_BIN *pf)
{
  INT i,j,min_key,min_index,finish=0;
  INT *buf_INT;
  FLOAT *buf_FLOAT;
  INT nints, nfloats;

  /* copy local state to array */
  pf->state[degree] = pf->local;

  /* update downtree states */
  for (i=0; i<degree; i++)
    if (!pf->valid_state[i]) {
      GetConcentrate(i,pf->state+i,sizeof(NS_DIM_PREFIX PFILE_STATE_BIN));
      pf->valid_state[i] = 1;
#ifdef LOCAL_DEBUG
      UserWriteF("receiving state from %d: nints=%d, nfloats=%d fst=%d, lst=%d\n",
                 i,pf->state[i].nints,pf->state[i].nfloats,pf->state[i].first_key,pf->state[i].last_key);
#endif
    }
  /* now we have a valid state from all downtree nodes	*/
  /* and ourself.											*/

  /* determine smallest key */
  min_key = PFILE_MAX_INT;
  for (j=0; j<=degree; j++) {
    i = (pf->robin+j) % (degree+1);             /* makes load balancing in case of equal keys */
    if ( ((pf->state[i].nints>0) || (pf->state[i].nfloats>0)) && (pf->state[i].first_key<min_key) ) {
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
  if (min_key == PFILE_MAX_INT) {
    if (me == master)
      finish = 1;                   /* now all are finished */
    else {
      Concentrate(&pf->local,sizeof(NS_DIM_PREFIX PFILE_STATE_BIN));                   /* has (nints && nfloats)==0 !  */
      GetSpread(&finish,sizeof(INT));                                            /* the global state */
    }

    /* now we have finish from uptree */
    if (finish)     {
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
  if (min_index==degree) {
    /* it is our own data */
#ifdef LOCAL_DEBUG
    UserWriteF("Writing my own buffer, nints=%d, nfloats=%d\n",pf->local.nints,pf->local.nfloats);
#endif
    buf_INT = pf->buf_INT;
    buf_FLOAT = pf->buf_FLOAT;
    nints = pf->local.nints;
    nfloats = pf->local.nfloats;
    pf->local.nints = pf->local.nfloats = 0;              /* buffer is empty again */
    pf->local.first_key = pf->local.last_key = 0;
  } else {
    /* fetch downtree data, we have only state ! */
    GetConcentrate(min_index,pf->buf_INT2,pf->state[min_index].nints*sizeof(INT));
    GetConcentrate(min_index,pf->buf_FLOAT2,pf->state[min_index].nfloats*sizeof(FLOAT));
#ifdef LOCAL_DEBUG
    UserWriteF("Writing buffer from %d, nints=%d, nfloats=%d\n",/
               min_index,pf->state[min_index].nints,pf->state[min_index].nfloats);
#endif
    buf_INT = pf->buf_INT2;
    buf_FLOAT = pf->buf_FLOAT2;
    nints = pf->state[min_index].nints;
    nfloats = pf->state[min_index].nfloats;
    pf->valid_state[min_index] = 0;             /* get new state next time */
  }
  /* now process data */
  if (me==master) {
    if (nints>0) fwrite(buf_INT, sizeof(INT), nints, pf->stream);
    if (nfloats>0) fwrite(buf_FLOAT, sizeof(FLOAT), nfloats, pf->stream);
  }
  else {
    Concentrate(pf->state+min_index,sizeof(NS_DIM_PREFIX PFILE_STATE_BIN));             /* send state */
    Concentrate(buf_INT,nints*sizeof(INT));                               /* send INT data  */
    Concentrate(buf_FLOAT,nfloats*sizeof(FLOAT));                         /* send FLOAT data  */
  }

  /* Note: While sending data uptree we don't expect to get	*/
  /* an answer from uptree!									*/
  return(0);       /* we are not finished yet */
}
#endif

#ifdef ModelP
static INT append_buffer_bin_INT (NS_DIM_PREFIX PFILE_BIN *pf, INT *values, int n, INT key)
{
  int i;
#ifdef LOCAL_DEBUG
  UserWriteF("appending %d INTs with key %d, nints=%d, lst=%d\n",n,key,
             pf->local.nints,pf->local.last_key);
#endif
  if (n+pf->local.nints<PFILE_BUFFER_SIZE)
  {
    /* initialize keys if buffer is empty */
    if (pf->local.nints == 0) pf->local.first_key = key;

    /* append to buffer */
    for (i=0; i<n; i++)
      pf->buf_INT[pf->local.nints+i]=values[i];
    pf->local.nints += n;

    /* update last key */
    pf->local.last_key = key;

    return(1);             /* append successful */
  }

  if (n>PFILE_BUFFER_SIZE)
  {
    PrintErrorMessage('E',"pfile_puts","string larger than buffer");
    return(1);             /* treat as appended ! */
  }

  return(0);       /* could not append */
}
#endif

/****************************************************************************/
/** \brief Write tagged sequence of integers to parallel binary file

   \param pf - pointer to parallel binary file
   \param values - integers to be written
   \param n - number of integers
   \param key - key for this sequence of integers

   Writes a tagged item to the output file. Tags should be globally
   unique and locally increasing within each segment. If tags
   are out of order, a warning is issued and the item is treated
   as in order.

   Any number of calls to pfile_tagged_write_INT may be issued by
   a single process.

   \return <ul>
   <li> 0 if ok </li>
   <li> >0 if any errors are encountered, this is a LOCAL state </li>
   </ul>

   \sa
   'pfile_sync_bin', 'pfile_open_bin', 'pfile_close_bin'
 */
/****************************************************************************/

INT NS_DIM_PREFIX pfile_tagged_write_INT (PFILE_BIN *pf, INT *values, int n, INT key)
{
#ifndef ModelP
  fwrite(values, sizeof(INT), n, (FILE *) pf);
#else
  /* check order of keys */
  if ( (pf->local.nints > 0) && (key < pf->local.last_key) )
  {
    PrintErrorMessage('W',"pfile_tagged_write_INT","keys locally not ordered");
    key = pf->local.last_key+1;             /* treat as in order ! */
  }

  if ( (pf->local.nints > 0) && (key > pf->local.last_key+1) )       /* not consecutive case */
    if (flush_buffer_bin(pf))
    {
      PrintErrorMessage('E',"pfile_tagged_write_INT","unxepected finish");
      return(1);
    }

  if (!append_buffer_bin_INT(pf,values,n,key))
  {
    if (flush_buffer_bin(pf))
    {
      PrintErrorMessage('E',"pfile_tagged_write_INT","unxepected finish");
      return(1);
    }
    append_buffer_bin_INT(pf,values,n,key);
  }
#endif
  return(0);
}

/****************************************************************************/
#ifdef ModelP
static INT append_buffer_bin_FLOAT (NS_DIM_PREFIX PFILE_BIN *pf, FLOAT *values, int n, INT key)
{
  int i;
#ifdef LOCAL_DEBUG
  UserWriteF("appending %d FLOATs with key %d, nfloats=%d, lst=%d\n",n,key,
             pf->local.nfloats,pf->local.last_key);
#endif
  if (n+pf->local.nfloats<PFILE_BUFFER_SIZE)
  {
    /* initialize keys if buffer is empty */
    if (pf->local.nfloats == 0) pf->local.first_key = key;

    /* append to buffer */
    for (i=0; i<n; i++)
      pf->buf_FLOAT[pf->local.nfloats+i]=values[i];
    pf->local.nfloats += n;

    /* update last key */
    pf->local.last_key = key;

    return(1);             /* append successful */
  }

  if (n>PFILE_BUFFER_SIZE)
  {
    PrintErrorMessage('E',"pfile_puts","string larger than buffer");
    return(1);             /* treat as appended ! */
  }

  return(0);       /* could not append */
}
#endif

INT NS_DIM_PREFIX pfile_tagged_write_FLOAT (PFILE_BIN *pf, FLOAT *values, int n, INT key)
{
#ifndef ModelP
  fwrite(values, sizeof(FLOAT), n, (FILE *) pf);
#else
  /* check order of keys */
  if ( (pf->local.nfloats > 0) && (key < pf->local.last_key) )
  {
    PrintErrorMessage('W',"pfile_tagged_write_INT","keys locally not ordered");
    key = pf->local.last_key+1;             /* treat as in order ! */
  }

  if ( (pf->local.nfloats > 0) && (key > pf->local.last_key+1) )       /* not consecutive case */
    if (flush_buffer_bin(pf))
    {
      PrintErrorMessage('E',"pfile_tagged_write_INT","unxepected finish");
      return(1);
    }

  if (!append_buffer_bin_FLOAT(pf,values,n,key))
  {
    if (flush_buffer_bin(pf))
    {
      PrintErrorMessage('E',"pfile_tagged_write_INT","unxepected finish");
      return(1);
    }
    append_buffer_bin_FLOAT(pf,values,n,key);
  }
#endif
  return(0);
}

/****************************************************************************/
#ifdef ModelP
static INT append_buffer_bin_BYTE (NS_DIM_PREFIX PFILE_BIN *pf, unsigned char *values, int n, INT key)
{
  int i;
#ifdef LOCAL_DEBUG
  UserWriteF("appending %d BYTEs with key %d, nbytes=%d, lst=%d\n",n,key,
             pf->local.nbytes,pf->local.last_key);
#endif
  if (n+pf->local.nbytes<PFILE_BUFFER_SIZE)
  {
    /* initialize keys if buffer is empty */
    if (pf->local.nbytes == 0) pf->local.first_key = key;

    /* append to buffer */
    for (i=0; i<n; i++)
      pf->buf_BYTE[pf->local.nbytes+i]=values[i];
    pf->local.nbytes += n;

    /* update last key */
    pf->local.last_key = key;

    return(1);             /* append successful */
  }

  if (n>PFILE_BUFFER_SIZE)
  {
    PrintErrorMessage('E',"pfile_puts","string larger than buffer");
    return(1);             /* treat as appended ! */
  }

  return(0);       /* could not append */
}
#endif

INT NS_DIM_PREFIX pfile_tagged_write_BYTE (PFILE_BIN *pf, unsigned char *values, int n, INT key)
{
#ifndef ModelP
  fwrite(values, sizeof(unsigned char), n, (FILE *) pf);
#else
  /* check order of keys */
  if ( (pf->local.nbytes > 0) && (key < pf->local.last_key) )
  {
    PrintErrorMessage('W',"pfile_tagged_write_INT","keys locally not ordered");
    key = pf->local.last_key+1;             /* treat as in order ! */
  }

  if ( (pf->local.nbytes > 0) && (key > pf->local.last_key+1) )       /* not consecutive case */
    if (flush_buffer_bin(pf))
    {
      PrintErrorMessage('E',"pfile_tagged_write_INT","unxepected finish");
      return(1);
    }

  if (!append_buffer_bin_BYTE(pf,values,n,key))
  {
    if (flush_buffer_bin(pf))
    {
      PrintErrorMessage('E',"pfile_tagged_write_INT","unxepected finish");
      return(1);
    }
    append_buffer_bin_BYTE(pf,values,n,key);
  }
#endif
  return(0);
}

/****************************************************************************/
/** \brief Indicate end of binary file segment

   \param pf - pointer to parallel binary file

   At the end of a segment this function has to be called by
   all processes. The function returns when all process
   have reached the end the current segment. Then the segment
   counter is increased.

   \return <ul>
   <li> 0 if ok </li>
   <li> >0 if any errors are encountered, this is a LOCAL state </li>
   </ul>

   \sa
   'pfile_open_bin', 'pfile_close_bin'
 */
/****************************************************************************/

INT NS_DIM_PREFIX pfile_sync_bin (PFILE_BIN *pf)
{
#ifdef ModelP
  while (!flush_buffer_bin(pf)) ;
#endif
  return(0);
}

/****************************************************************************/
/** \brief Indicate end of binary file

   \param pf - pointer to parallel binary file

   At the end of the file this function has to be called by
   all processes. The function returns when all process
   have reached the end of file. Then the file is closed
   and buffer space is released.

   \return <ul>
   <li> 0 if ok </li>
   <li> >0 if any errors are encountered, this is a LOCAL state </li>
   </ul>

   \sa
   'pfile_open_bin', 'pfile_sync_bin'
 */
/****************************************************************************/

INT NS_DIM_PREFIX pfile_close_bin (PFILE_BIN *pf)
{
#ifndef ModelP
  fclose( (FILE *) pf);
#else
  /* wait until all processors reach end of segment */
  while (!flush_buffer_bin(pf)) ;
  /* close output file */
  if (me == master) fclose(pf->stream);
  /* free data structure */
  free(pf);
#endif
  return(0);
}
