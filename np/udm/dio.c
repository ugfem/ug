// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  dio.c															*/
/*																			*/
/* Purpose:   input/output of data			                                                                */
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   16.12.96 begin,												*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdio.h>
#include "bio.h"
#include "dio.h"

#define __MGIO_USE_IN_UG__

#ifdef __MGIO_USE_IN_UG__

        #include "defaults.h"
        #include "fileopen.h"
        #include "np.h"

#endif

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define DIO_TITLE_LINE                          "####.sparse.data.storage.format.####"

#define DIO_INTSIZE                                     100
#define DIO_BUFFERSIZE                          128

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



/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static FILE *stream;                            /* file                                                 */
static int datapaths_set;                       /* pathes used in ug			*/
static char buffer[DIO_BUFFERSIZE]; /* general purpose buffer		*/
static int intList[DIO_INTSIZE];        /* general purpose integer list */

/* RCS string */
#ifdef __MGIO_USE_IN_UG__

RCSID("$Header$",UG_RCS_STRING)

#endif

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   Read_OpenDTFile - opens file for reading data

   SYNOPSIS:
   int Read_OpenDTFile (char *filename);

   PARAMETERS:
   .  filename - name of file

   DESCRIPTION:
   opens a file with specified name for reading

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 when error occured.

   SEE ALSO:
   D*/
/****************************************************************************/

int Read_OpenDTFile (char *filename)
{

#ifdef __MGIO_USE_IN_UG__
  if (datapaths_set) stream = FileOpenUsingSearchPaths(filename,"r","datapaths");
  else stream = fileopen(filename,"r");
#else
  stream = fopen(filename,"r");
#endif

  if (stream==NULL) return (1);

  return (0);
}

/****************************************************************************/
/*D
   Write_OpenDTFile - opens file for reading

   SYNOPSIS:
   int Write_OpenDTFile (char *filename);

   PARAMETERS:
   .  filename - name of file

   DESCRIPTION:
   opens a file with specified name for writing data

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 when error occured.

   SEE ALSO:
   D*/
/****************************************************************************/

int Write_OpenDTFile (char *filename)
{

#ifdef __MGIO_USE_IN_UG__
  if (datapaths_set) stream = FileOpenUsingSearchPaths(filename,"w","datapaths");
  else stream = fileopen(filename,"w");
#else
  stream = fopen(filename,"w");
#endif

  if (stream==NULL) return (1);
  return (0);
}

/****************************************************************************/
/*
   Read_DT_General - reads general information about mg

   SYNOPSIS:
   int Read_MG_General (MGIO_MG_GENERAL *mg_general);

   PARAMETERS:
   .  mg_general - general information about mg

   DESCRIPTION:
   function reads general information about the mg

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Read_DT_General (DIO_GENERAL *dio_general)
{
  int i;

  /* initialize basic i/o */
  if (Bio_Initialize(stream,BIO_ASCII)) return (1);

  /* head always in ACSII */
  if (Bio_Read_string(buffer)) return (1);if (strcmp(buffer,DIO_TITLE_LINE)!=0) return (1);
  if (Bio_Read_mint(1,intList)) return (1);
  dio_general->mode               = intList[0];

  /* re-initialize basic i/o */
  if (Bio_Initialize(stream,dio_general->mode)) return (1);

  /* now special mode */
  if (Bio_Read_string(dio_general->version)) return (1);
  if (Bio_Read_mint(2,intList)) return (1);
  dio_general->magic_cookie       = intList[0];
  dio_general->nVD                        = intList[1];
  for (i=0; i<dio_general->nVD; i++)
  {
    if (Bio_Read_string(dio_general->VDname[i])) return (1);
    if (Bio_Read_mint(1,dio_general->VDncomp+i)) return (1);
  }
  if (Bio_Read_mint(1,intList)) return (1);
  dio_general->ndata              = intList[0];

  return (0);
}

/****************************************************************************/
/*
   Write_MG_General - writes general information about mg

   SYNOPSIS:
   int Write_MG_General (MGIO_MG_GENERAL *mg_general);

   PARAMETERS:
   .  mg_general - general information about mg

   DESCRIPTION:
   function writes general information about the mg

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int     Write_DT_General (DIO_GENERAL *dio_general)
{
  int i;

  /* initialize basic i/o */
  if (Bio_Initialize(stream,BIO_ASCII)) return (1);

  /* head always in ACSII */
  if (Bio_Write_string(DIO_TITLE_LINE)) return (1);
  intList[0] = dio_general->mode;
  if (Bio_Write_mint(1,intList)) return (1);

  /* re-initialize basic i/o */
  if (Bio_Initialize(stream,dio_general->mode)) return (1);

  /* now special mode */
  if (Bio_Write_string(dio_general->version)) return (1);
  intList[0] = dio_general->magic_cookie;
  intList[1] = dio_general->nVD;
  if (Bio_Write_mint(2,intList)) return (1);
  for (i=0; i<dio_general->nVD; i++)
  {
    if (Bio_Write_string(dio_general->VDname[i])) return (1);
    if (Bio_Write_mint(1,dio_general->VDncomp+i)) return (1);
  }
  intList[0] = dio_general->ndata;
  if (Bio_Write_mint(1,intList)) return (1);

  return (0);
}

/****************************************************************************/
/*
   CloseFile - close the file

   SYNOPSIS:
   int CloseFile ();

   PARAMETERS:

   DESCRIPTION:
   close the file

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int CloseDTFile ()
{
  if (fclose(stream)!=NULL) return (1);
  return (0);
}

/****************************************************************************/
/*
   DIO_Init - init input/output for data

   SYNOPSIS:
   int DIO_Init (void);

   PARAMETERS:

   DESCRIPTION:
   init the i/o of data

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if read error.

   SEE ALSO:
 */
/****************************************************************************/

int DIO_Init (void)
{

#ifdef __MGIO_USE_IN_UG__

  /* path to grid-dirs */
  datapaths_set = 0;
  if (ReadSearchingPaths(DEFAULTSFILENAME,"datapaths")==0)
    datapaths_set = 1;

#endif

  return (0);
}
