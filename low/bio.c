// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*! \file bio.c
 * \ingroup low
 */

/** \addtogroup low
 *
 * @{
 */

/****************************************************************************/
/*                                                                          */
/* File:      bio.c                                                         */
/*                                                                          */
/* Purpose:   basic input/output                                            */
/*                                                                          */
/* Author:    Klaus Johannsen                                               */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70550 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
/*                                                                          */
/* History:   09.12.96 begin,                                               */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/* system include files                                                     */
/* application include files                                                */
/*                                                                          */
/****************************************************************************/

#include <config.h>
#include <stdio.h>
#include <string.h>

#include <rpc/rpc.h>    /* to include xdr.h in a portable way */

#include "general.h"

#include "bio.h"

USING_UG_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*    compile time constants defining static data size (i.e. arrays)        */
/*    other constants                                                       */
/*    macros                                                                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/* in the corresponding include file!)                                      */
/*                                                                          */
/****************************************************************************/

typedef int (*R_mint_proc)(int n, int *intList);
typedef int (*W_mint_proc)(int n, const int *intList);
typedef int (*R_mdouble_proc)(int n, double *doubleList);
typedef int (*W_mdouble_proc)(int n, const double *doubleList);
typedef int (*R_string_proc)(char *string);
typedef int (*W_string_proc)(const char *string);

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* file */
static FILE *stream;
static int n_byte;
static fpos_t pos;
static XDR xdrs;

/* low level read/write functions */
static R_mint_proc Read_mint;
static W_mint_proc Write_mint;
static R_mdouble_proc Read_mdouble;
static W_mdouble_proc Write_mdouble;
static R_string_proc Read_string;
static W_string_proc Write_string;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* debug i/o                                                                */
/*                                                                          */
/****************************************************************************/

static int XDR_Read_mint (int n, int *intList)
{
  int i;

  for (i=0; i<n; i++)
  {
    /*if (!xdr_int(&xdrs,&(intList[i]))) return (1); quick hack to cure Cray's xdr_int error Christian Wrobel 980529*/
    if (!xdr_u_int(&xdrs,(unsigned int*)&(intList[i]))) return (1);
#if defined __T3D__ || defined __T3E__
    intList[i] = (short int) intList[i];            /* sign extension due to xdr_u_int, 4/8 Byte conversion and bug in Cray's xdr on certain platforms; Christian Wrobel 980529 */
#endif
  }
  return (0);
}

static int XDR_Write_mint (int n, const int *intList)
{
  int i;

  for (i=0; i<n; i++)
  {
    /*if (!xdr_int(&xdrs,&(intList[i]))) return (1); quick hack to cure Cray's
       xdr_int error Christian Wrobel 980529*/
    if (!xdr_u_int(&xdrs,(unsigned int*)&(intList[i]))) return (1);
    n_byte += 4;
  }
  return (0);
}

static int XDR_Read_mdouble (int n, double *doubleList)
{
  int i;

  for (i=0; i<n; i++)
    if (!xdr_double(&xdrs,&(doubleList[i]))) return (1);
  return (0);
}

static int XDR_Write_mdouble (int n, const double *doubleList)
{
  int i;

  for (i=0; i<n; i++)
  {
    if (!xdr_double(&xdrs,(double*) &(doubleList[i]))) return (1);
    n_byte += 8;
  }
  return (0);
}

static int XDR_Read_string (char *string)
{
  int i,len;

  if (fscanf(stream,"%d ",&len)!=1) return (1);
  for (i=0; i<len; i++)
  {
    string[i] = fgetc(stream);
    if (string[i]==EOF)
      return (1);
  }
  string[i] = fgetc(stream);
  if (string[i]!='\n') return (1);
  string[i] = '\0';

  return (0);
}

static int XDR_Write_string (const char *string)
{
  int i,m,len;

  len = strlen(string);
  m = fprintf(stream,"%d ",len);
  if (m<0) return (1);
  n_byte += m;
  for (i=0; i<len; i++)
    if (fputc(string[i],stream)==EOF)
      return (1);
  m = fprintf(stream,"\n");
  if (m<0) return (1);
  n_byte += len+m;

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* ascii i/o                                                                */
/*                                                                          */
/****************************************************************************/


static int ASCII_Read_mint (int n, int *intList)
{
  int i;

  for (i=0; i<n; i++)
    if (fscanf(stream,"%d\n",intList+i)!=1) return (1);
  return (0);
}

static int ASCII_Write_mint (int n, const int *intList)
{
  int i,m;

  for (i=0; i<n; i++)
  {
    m = fprintf(stream,"%d\n",intList[i]);
    if (m<0) return (1);
    n_byte += m;
  }
  return (0);
}

static int ASCII_Read_mdouble (int n, double *doubleList)
{
  int i;

  for (i=0; i<n; i++)
    if (fscanf(stream,"%lg\n",doubleList+i)!=1) return (1);
  return (0);
}

static int ASCII_Write_mdouble (int n, const double *doubleList)
{
  int i,m;

  for (i=0; i<n; i++)
  {
    m = fprintf(stream,"%g\n",doubleList[i]);
    if (m<0) return (1);
    n_byte += m;
  }
  return (0);
}

static int ASCII_Read_string (char *string)
{
  int i,len;

  if (fscanf(stream,"%d\n",&len)!=1) return (1);
  for (i=0; i<len; i++)
  {
    string[i] = fgetc(stream);
    if (string[i]==EOF)
      return (1);
  }
  string[i] = fgetc(stream);
  if (string[i]!=' ') return (1);
  string[i] = '\0';

  return (0);
}


static int ASCII_Write_string (const char *string)
{
  int i,m,len;

  len = strlen(string);
  m = fprintf(stream,"%d\n",len);
  if (m<0) return (1);
  n_byte += m;
  for (i=0; i<len; i++)
    if (fputc(string[i],stream)==EOF)
      return (1);
  m = fprintf(stream," ");
  if (m<0) return (1);
  n_byte += len+m;

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* binary i/o                                                               */
/*                                                                          */
/****************************************************************************/


static int BIN_Read_mint (int n, int *intList)
{
  if (fread((void*)intList,sizeof(int)*n,1,stream)!=1) return (1);
  return (0);
}

static int BIN_Write_mint (int n, const int *intList)
{
  if (fwrite((void*)intList,sizeof(int)*n,1,stream)!=1) return (1);
  n_byte += n*sizeof(int);
  return (0);
}

static int BIN_Read_mdouble (int n, double *doubleList)
{
  if (fread((void*)doubleList,sizeof(double)*n,1,stream)!=1) return (1);
  return (0);
}

static int BIN_Write_mdouble (int n, const double *doubleList)
{
  if (fwrite((void*)doubleList,sizeof(double)*n,1,stream)!=1) return (1);
  n_byte += n*sizeof(double);
  return (0);
}

static int BIN_Read_string (char *string)
{
  int i,len;

  if (fscanf(stream,"%d ",&len)!=1) return (1);
  for (i=0; i<len; i++)
  {
    string[i] = fgetc(stream);
    if (string[i]==EOF)
      return (1);
  }
  string[i] = fgetc(stream);
  if (string[i]!=' ') return (1);
  string[i] = '\0';

  return (0);
}


static int BIN_Write_string (const char *string)
{
  int i,m,len;

  len = strlen(string);
  m = fprintf(stream,"%d ",len);
  if (m<0) return (1);
  n_byte += m;
  for (i=0; i<len; i++)
    if (fputc(string[i],stream)==EOF)
      return (1);
  m = fprintf(stream," ");
  if (m<0) return (1);
  n_byte += len+m;

  return (0);
}
/****************************************************************************/
/*                                                                          */
/* exported i/o                                                             */
/*                                                                          */
/****************************************************************************/

int NS_PREFIX Bio_Initialize (FILE *file, int mode, char rw)
{
  stream = file;

  switch (mode)
  {
  case BIO_XDR :
    if (rw=='r') xdrstdio_create(&xdrs,file,XDR_DECODE);
    else if (rw=='w') xdrstdio_create(&xdrs,file,XDR_ENCODE);
    else return (1);
    Read_mint       = XDR_Read_mint;
    Read_mdouble = XDR_Read_mdouble;
    Read_string = XDR_Read_string;
    Write_mint      = XDR_Write_mint;
    Write_mdouble = XDR_Write_mdouble;
    Write_string = XDR_Write_string;
    break;
  case BIO_ASCII :
    Read_mint       = ASCII_Read_mint;
    Read_mdouble = ASCII_Read_mdouble;
    Read_string = ASCII_Read_string;
    Write_mint      = ASCII_Write_mint;
    Write_mdouble = ASCII_Write_mdouble;
    Write_string = ASCII_Write_string;
    break;
  case BIO_BIN :
    Read_mint       = BIN_Read_mint;
    Read_mdouble = BIN_Read_mdouble;
    Read_string = BIN_Read_string;
    Write_mint      = BIN_Write_mint;
    Write_mdouble = BIN_Write_mdouble;
    Write_string = BIN_Write_string;
    break;
  default :
    return (1);
  }

  return (0);
}

int NS_PREFIX Bio_Read_mint (int n, int *intList)
{
  return ((*Read_mint)(n,intList));
}

int NS_PREFIX Bio_Write_mint (int n, int *intList)
{
  return ((*Write_mint)(n,intList));
}

int NS_PREFIX Bio_Read_mdouble (int n, double *doubleList)
{
  return ((*Read_mdouble)(n,doubleList));
}

int NS_PREFIX Bio_Write_mdouble (int n, double *doubleList)
{
  return ((*Write_mdouble)(n,doubleList));
}

int NS_PREFIX Bio_Read_string (char *string)
{
  return ((*Read_string)(string));
}

int NS_PREFIX Bio_Write_string (const char *string)
{
  return ((*Write_string)(string));
}

int NS_PREFIX Bio_Jump_From (void)
{
  n_byte = 0;
  if (fgetpos(stream,&pos)) return (1);
  if (fprintf(stream," %20d ",n_byte)<0) return (1);

  return (0);
}

int NS_PREFIX Bio_Jump_To (void)
{
  fpos_t act;

  if (fgetpos(stream,&act)) return (1);
  if (fsetpos(stream,&pos)) return (1);
  if (fprintf(stream," %20d ",n_byte)<0) return (1);
  if (fsetpos(stream,&act)) return (1);

  return (0);
}

int NS_PREFIX Bio_Jump (int dojump)
{
  int jump;

  if (fscanf(stream," %20d ",&jump)!=1) return (1);
  if (dojump==0) return (0);
  while(jump>0)
  {
    if (fgetc(stream)==EOF) return (1);
    jump--;
  }

  return (0);
}

/** @} */
