// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  scan.c	                                                                                                */
/*																			*/
/* Purpose:   tools for reading script arguments                                */
/*																			*/
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   November 23, 1996                                                                         */
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "general.h"
#include "debug.h"

#include "gm.h"
#include "ugenv.h"
#include "devices.h"

#include "formats.h"
#include "pcr.h"
#include "numproc.h"
#include "np.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define OPTIONLEN                       32
#define OPTIONLENSTR            "31"
#define VALUELEN                        64
#define VALUELENSTR                     "63"

/* token seperators for ReadVecType... */
#define TYPESEP                 "|"
#define COMPSEP                 " \t:"

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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   ReadArgvDOUBLE - Read command strings

   SYNOPSIS:
   INT ReadArgvDOUBLE (const char *name, DOUBLE *a, INT argc, char **argv);

   PARAMETERS:
   .  name - name of the argument
   .  a - DOUBLE value
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads command strings and returns a DOUBLE value.

   RETURN VALUE:
   INT
   .n    0 if the argument was found and a DOUBLE value could be read
   .n    1 else.
   D*/
/****************************************************************************/

INT ReadArgvDOUBLE (const char *name, DOUBLE *a, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN];
  double value;

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if (sscanf(argv[i],"%s %lf",option,&value)!=2)
        continue;
      if (strcmp(option,name) == 0)
      {
        a[0] = value;
        return(0);
      }
    }

  REP_ERR_RETURN (1);
}

/****************************************************************************/
/*D
   ReadArgvINT - Read command strings

   SYNOPSIS:
   INT ReadArgvINT (const char *name, INT *j, INT argc, char **argv);

   PARAMETERS:
   .  name - name of the argument
   .  j - integer value
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads command strings and returns an integer value.

   RETURN VALUE:
   INT
   .n    0 if the argument was found and a value could be read
   .n    1 else.
   D*/
/****************************************************************************/

INT ReadArgvINT (const char *name, INT *j, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN];
  int value;

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if (sscanf(argv[i],"%s %d",option,&value)!=2)
        continue;
      if (strcmp(option,name) == 0)
      {
        j[0] = value;
        return(0);
      }
    }

  REP_ERR_RETURN (1);
}

/****************************************************************************/
/*D
   ReadArgvChar - Read command strings

   SYNOPSIS:
   INT ReadArgvChar (const char *name, char *buffer, INT argc, char **argv);

   PARAMETERS:
   .  name - name of the argument
   .  buffer - string
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads command strings and returns an integer value.

   RETURN VALUE:
   INT
   .n    0 if the argument was found and a value could be read
   .n    1 else.
   D*/
/****************************************************************************/

INT ReadArgvChar (const char *name, char *buffer, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN];
  char value[VALUELEN];

  strcpy(buffer,"");
  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0]) {
      if (sscanf(argv[i],
                 expandfmt(CONCAT5("%",OPTIONLENSTR,"[a-zA-Z0-9_] %",
                                   VALUELENSTR,"[ -~]")),option,value)!=2)
        continue;
      if (strcmp(option,name) == 0) {
        strcpy(buffer,value);
        return(0);
      }
    }

  REP_ERR_RETURN (1);
}

/****************************************************************************/
/*D
   ReadArgvDisplay - Read command strings

   SYNOPSIS:
   INT ReadArgvDisplay (INT argc, char **argv);

   PARAMETERS:
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads the display status.

   RETURN VALUE:
   INT
   .n    PCR_NO_DISPLAY     no display (default if not specified)
   .n    PCR_RED_DISPLAY    reduced display
   .n    PCR_FULL_DISPLAY   full display
   D*/
/****************************************************************************/

INT ReadArgvDisplay (INT argc, char **argv)
{
  INT i;
  char value[VALUELEN];

  for (i=0; i<argc; i++)
    if (strncmp(argv[i],"display",7)==0)
    {
      if (sscanf(argv[i],"display %s",value) != 1)
        continue;
      if (strcmp(value,"no") == 0)
        return(PCR_NO_DISPLAY);
      else if (strcmp(value,"red") == 0)
        return(PCR_RED_DISPLAY);
      else if (strcmp(value,"full") == 0)
        return(PCR_FULL_DISPLAY);
    }

  return(PCR_NO_DISPLAY);
}

/****************************************************************************/
/*D
   ReadArgvOption - Read command strings

   SYNOPSIS:
   INT ReadArgvOption (const char *name, INT argc, char **argv);

   PARAMETERS:
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads command strings and returns an integer value.

   RETURN VALUE:
   INT
   .n    0 if the option is not set
   .n    1 else.
   D*/
/****************************************************************************/

INT ReadArgvOption (const char *name, INT argc, char **argv)
{
  INT i;
  char option[OPTIONLEN];
  int value;

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if (sscanf(argv[i],"%s %d",option,&value) == 2)
        if (strcmp(option,name) == 0)
          return(value);
      if (strcmp(argv[i],name) == 0)
        return(1);
    }

  return (0);
}

/****************************************************************************/
/*D
   ReadArgvPosition - Read command strings

   SYNOPSIS:
   INT ReadArgvPosition (const char *name, INT argc, char **argv, DOUBLE *pos);

   PARAMETERS:
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector
   .  pos - position vector

   DESCRIPTION:
   This function reads command strings and returns a position.

   RETURN VALUE:
   INT
   .n    0 if the argument was found and a position could be read
   .n    1 else.
   D*/
/****************************************************************************/

INT ReadArgvPosition (const char *name, INT argc, char **argv, DOUBLE *pos)
{
  INT i;
  char option[OPTIONLEN];
  float x,y,z;

  for (i=0; i<argc; i++)
    if (argv[i][0]==name[0])
    {
      if (sscanf(argv[i],"%s %f %f %f",option,&x,&y,&z)!=DIM+1)
        continue;
      if (strcmp(option,name) == 0)
      {
        pos[0] = x;
        pos[1] = y;
              #ifdef __THREEDIM__
        pos[2] = z;
                      #endif
        return(0);
      }
    }

  REP_ERR_RETURN(1);
}

/****************************************************************************/
/*D
   ReadArgvVecDesc - Read command strings

   SYNOPSIS:
   VECDATA_DESC *ReadArgvVecDesc (MULTIGRID *theMG, const char *name,
                                  INT argc, char **argv);

   PARAMETERS:
   .  theMG - pointer to a multigrid
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads a symbol name from the command strings and returns
   a pointer to the corresponding vector descriptor.

   CAUTION: If no template is specified the first vector template is used.

   This call locks the vector descriptor for dynamic allocation.

   SYNTAX:
   $s <vec desc name>[/<template name>]

   RETURN VALUE:
   VECDATA_DESC *
   .n    pointer to vector descriptor
   .n    NULL if error occurs
   D*/
/****************************************************************************/

VECDATA_DESC *ReadArgvVecDesc (MULTIGRID *theMG, const char *name,
                               INT argc, char **argv)
{
  VECDATA_DESC *vd;
  char value[VALUELEN],vdname[NAMESIZE],tname[NAMESIZE];
  INT res;

  if (ReadArgvChar(name,value,argc,argv))
    REP_ERR_RETURN (NULL);

  res = sscanf(value,expandfmt(CONCAT5("%",NAMELENSTR,"[a-zA-Z0-9_] / %",NAMELENSTR,"[a-zA-Z0-9_]")),vdname,tname);
  vd = GetVecDataDescByName(theMG,vdname);
  if (vd == NULL)
  {
    if (res==2)
      vd = CreateVecDescOfTemplate (theMG,vdname,tname);
    else
      /* taking default template */
      vd = CreateVecDescOfTemplate (theMG,vdname,NULL);
  }
  if (vd == NULL) REP_ERR_RETURN (NULL);

  VM_LOCKED(vd) = 1;

  return(vd);
}

/****************************************************************************/
/*D
   ReadArgvMatDesc - Read command strings

   SYNOPSIS:
   MATDATA_DESC *ReadArgvMatDesc (MULTIGRID *theMG, const char *name,
   INT argc, char **argv);

   PARAMETERS:
   .  theMG - pointer to a multigrid
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads a symbol name from the command strings and returns
   a pointer to the corresponding matrix descriptor.
   This call locks the matrix descriptor for dynamic allocation.

   RETURN VALUE:
   MATDATA_DESC *
   .n    pointer to matrix descriptor
   .n    NULL if error occurs
   D*/
/****************************************************************************/

MATDATA_DESC *ReadArgvMatDesc (MULTIGRID *theMG, const char *name,
                               INT argc, char **argv)
{
  MATDATA_DESC *md;
  char value[VALUELEN],mdname[NAMESIZE],tname[NAMESIZE];
  INT res;

  if (ReadArgvChar(name,value,argc,argv))
    REP_ERR_RETURN (NULL);

  res = sscanf(value,"%s/%s",mdname,tname);
  md = GetMatDataDescByName(theMG,mdname);
  if (md == NULL)
  {
    if (res==2)
      md = CreateMatDescOfTemplate (theMG,mdname,tname);
    else
      /* taking default template */
      md = CreateMatDescOfTemplate (theMG,mdname,NULL);
  }
  if (md == NULL) REP_ERR_RETURN (NULL);

  VM_LOCKED(md) = 1;

  return(md);
}

/****************************************************************************/
/*D
   ReadArgvNumProc - Read command strings

   SYNOPSIS:
   NP_BASE *ReadArgvNumProc (MULTIGRID *theMG, const char *name, char *class,
   INT argc, char **argv);

   PARAMETERS:
   .  theMG - pointer to a multigrid
   .  name - name of the argument
   .  class - name of the class
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads a num proc name from the command strings and returns
   a pointer to the num proc.

   RETURN VALUE:
   NP_BASE *
   .n    pointer to num proc
   .n    NULL if error occurs
   D*/
/****************************************************************************/

NP_BASE *ReadArgvNumProc (MULTIGRID *theMG, const char *name, const char *class,
                          INT argc, char **argv)
{
  char value[VALUELEN];

  if (ReadArgvChar(name,value,argc,argv))
    REP_ERR_RETURN (NULL);

  return(GetNumProcByName(theMG,value,class));
}

/****************************************************************************/
/*
   ReadVecTypeINTs - Read a number of INTs from the input string

   SYNOPSIS:
   INT ReadVecTypeINTs (char *str, INT n, INT nINT[MAXVECTORS],
   INT theINTs[][MAXVECTORS]);

   PARAMETERS:
   .  str - input string
   .  n - maximal number of INTs
   .  nINT[MAXVECTORS] - number per vector type
   .  theINTs[][MAXVECTORS] - array to store the numbers

   DESCRIPTION:
   This function reads a number of integer values from the input string.
   It is used for the init routines of numprocs to transfor variables
   from the shell to the numproc. Then, 'str' is one argument of npinit
   of the format

   .  <sc~int~list>  - [nd <int  list>] | [ed <int list>] | [el <int list>] | [si <int list>]
   .  <int~list>  - [<int>[:<int>]*]

   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata


   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    else if error occured.
 */
/****************************************************************************/

INT ReadVecTypeINTs (char *str, INT n, INT nINT[MAXVECTORS], INT theINTs[][MAXVECTORS])
{
  char *s,*tok,*typetok[MAXVECTORS];
  INT type;
  int iValue;

  for (type=0; type<MAXVECTORS; type++) {nINT[type] = 0; typetok[type] = NULL;}

  for (tok=strtok(str,TYPESEP); tok!=NULL; tok=strtok(NULL,TYPESEP))
    if ((s=strstr(tok,"nd"))!=NULL)
      typetok[NODEVECTOR] = s+2;
    else if ((s=strstr(tok,"ed"))!=NULL)
      typetok[EDGEVECTOR] = s+2;
    else if ((s=strstr(tok,"el"))!=NULL)
      typetok[ELEMVECTOR] = s+2;
                #ifdef __THREEDIM__
  else if ((s=strstr(tok,"si"))!=NULL)
    typetok[SIDEVECTOR] = s+2;
                #endif
  else REP_ERR_RETURN (1);

  for (type=0; type<MAXVECTORS; type++)
    if (typetok[type]!=NULL)
      for (tok=strtok(typetok[type],COMPSEP); tok!=NULL; tok=strtok(NULL,COMPSEP))
      {
        if (nINT[type]>=n) REP_ERR_RETURN (2);

        if (sscanf(tok,"%d",&iValue)!=1)
          REP_ERR_RETURN (3)
          else
            theINTs[nINT[type]++][type] = (INT) iValue;
      }

  return (NUM_OK);
}

/****************************************************************************/
/*
   ReadVecTypeDOUBLEs - Read a number of DOUBLEs from the input string

   SYNOPSIS:
   INT ReadVecTypeDOUBLEs (char *str, INT n, INT nDOUBLE[MAXVECTORS],
   DOUBLE theDOUBLEs[][MAXVECTORS]);

   PARAMETERS:
   .  str - input string
   .  n - maximal number of DOUBLEs
   .  nDOUBLE[MAXVECTORS] - number per vector type
   .  theDOUBLEs[][MAXVECTORS] - array to store the numbers

   DESCRIPTION:
   This function reads a number of double values from the input string.
   It is used for the init routines of numprocs to transform variables
   from the shell to the numproc. Then, 'str' is one argument of npinit
   of the format

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - [<double>[:<double>]*]
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    else if error occured.
 */
/****************************************************************************/

INT ReadVecTypeDOUBLEs (char *str, INT n, INT nDOUBLE[MAXVECTORS], DOUBLE theDOUBLEs[][MAXVECTORS])
{
  char *s,*tok,*typetok[MAXVECTORS],*notypetok;
  INT type,found;
  double lfValue;

  for (type=0; type<MAXVECTORS; type++) {nDOUBLE[type] = 0; typetok[type] = NULL;}
  notypetok = NULL;

  for (tok=strtok(str,TYPESEP); tok!=NULL; tok=strtok(NULL,TYPESEP))
    if ((s=strstr(tok,"nd"))!=NULL)
      typetok[NODEVECTOR] = s+2;
    else if ((s=strstr(tok,"ed"))!=NULL)
      typetok[EDGEVECTOR] = s+2;
    else if ((s=strstr(tok,"el"))!=NULL)
      typetok[ELEMVECTOR] = s+2;
                #ifdef __THREEDIM__
  else if ((s=strstr(tok,"sd"))!=NULL)
    typetok[SIDEVECTOR] = s+2;
                #endif
  else
    notypetok = tok;

  found = 0;
  for (type=0; type<MAXVECTORS; type++)
    if (typetok[type]!=NULL)
      for (tok=strtok(typetok[type],COMPSEP); tok!=NULL; tok=strtok(NULL,COMPSEP))
      {
        found++;
        if (nDOUBLE[type]>=n) REP_ERR_RETURN (2);

        if (sscanf(tok,"%lf",&lfValue)!=1)
          REP_ERR_RETURN (3)
          else
            theDOUBLEs[nDOUBLE[type]++][type] = (DOUBLE) lfValue;
      }

  if (notypetok!=NULL)
  {
    if (found)
      REP_ERR_RETURN (NUM_ERROR);

    /* there is only one token witout type label */

    /* is there only one value? */
    found = 0;
    for (tok=strtok(notypetok,COMPSEP); tok!=NULL; tok=strtok(NULL,COMPSEP))
      found++;
    if (found!=1)
      REP_ERR_RETURN (NUM_ERROR)
      else
        REP_ERR_RETURN (NUM_TYPE_MISSING);
  }

  return (NUM_OK);
}

/****************************************************************************/
/*
   ReadVecTypeOrder - Read a number of INTs from the input string

   SYNOPSIS:
   INT ReadVecTypeOrder (char *str, INT n, INT MaxPerType,
   INT *nOrder, INT theOrder[]);

   PARAMETERS:
   .  str - input string
   .  n - maximal number
   .  MaxPerType - maximal number per type
   .  nOrder - number per type
   .  theOrder[] - array to store the integers

   DESCRIPTION:
   This function reads a number of INTs from the input string.
   It is used in an equation block smoother.

   SEE ALSO:
   ebgs

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    else if error occured.
 */
/****************************************************************************/

INT ReadVecTypeOrder (char *str, INT n, INT MaxPerType, INT *nOrder, INT theOrder[])
{
  char *token;
  INT ni;
  int iValue;

  ni = 0;
  for (token=strtok(str,COMPSEP); token!=NULL; token=strtok(NULL,COMPSEP))
  {
    if (ni>=n) REP_ERR_RETURN (1);

    if              ((sscanf(token,"nd%d",&iValue)==1) && (iValue<MaxPerType))
      theOrder[ni++] = NODEVECTOR*MaxPerType + (INT) iValue;
    else if ((sscanf(token,"ed%d",&iValue)==1) && (iValue<MaxPerType))
      theOrder[ni++] = EDGEVECTOR*MaxPerType + (INT) iValue;
    else if ((sscanf(token,"el%d",&iValue)==1) && (iValue<MaxPerType))
      theOrder[ni++] = ELEMVECTOR*MaxPerType + (INT) iValue;
                #ifdef __THREEDIM__
    else if ((sscanf(token,"si%d",&iValue)==1) && (iValue<MaxPerType))
      theOrder[ni++] = SIDEVECTOR*MaxPerType + (INT) iValue;
                #endif
    else REP_ERR_RETURN (2);
  }

  *nOrder = ni;

  return (NUM_OK);
}

/****************************************************************************/
/*
   ReadVecTypeNUMPROCs - Read a number of NUMPROCs from the input string

   SYNOPSIS:
   INT ReadVecTypeNUMPROCs (char *str, INT n, INT nNUMPROC[MAXVECTORS],
   NUM_PROC *theNUMPROCs[][MAXVECTORS]);

   PARAMETERS:
   .  str - input string
   .  n - maximal number of blocks
   .  nNUMPROC[MAXVECTORS] - number per vector type
   .  theNUMPROCs[][MAXVECTORS] - array to store the numprocs

   DESCRIPTION:
   This function reads a number of NUMPROCs from the input string.
   It is used to read the different smoothers per block in an
   equation block smoother.

   SEE ALSO:
   ebgs

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    else if error occured.
 */
/****************************************************************************/

INT ReadVecTypeNUMPROCs (MULTIGRID *theMG, char *str, char *class_name, INT n, INT nNUMPROC[MAXVECTORS], NP_BASE *theNUMPROCs[][MAXVECTORS])
{
  char *s,*tok,*typetok[MAXVECTORS];
  INT type;

  for (type=0; type<MAXVECTORS; type++) {nNUMPROC[type] = 0; typetok[type] = NULL;}

  for (tok=strtok(str,TYPESEP); tok!=NULL; tok=strtok(NULL,TYPESEP))
    if ((s=strstr(tok,"nd"))!=NULL)
      typetok[NODEVECTOR] = s+2;
    else if ((s=strstr(tok,"ed"))!=NULL)
      typetok[EDGEVECTOR] = s+2;
    else if ((s=strstr(tok,"el"))!=NULL)
      typetok[ELEMVECTOR] = s+2;
                #ifdef __THREEDIM__
  else if ((s=strstr(tok,"si"))!=NULL)
    typetok[SIDEVECTOR] = s+2;
                #endif
  else REP_ERR_RETURN (1);

  for (type=0; type<MAXVECTORS; type++)
    if (typetok[type]!=NULL)
      for (tok=strtok(typetok[type],COMPSEP); tok!=NULL; tok=strtok(NULL,COMPSEP))
      {
        if (nNUMPROC[type]>=n) REP_ERR_RETURN (2);

        if ((theNUMPROCs[nNUMPROC[type]++][type]=GetNumProcByName(theMG,tok,class_name))==NULL)
          REP_ERR_RETURN (3);
      }

  return (NUM_OK);
}

/****************************************************************************/
/*D
   sc_cmp - Compare VEC_SCALARs

   SYNOPSIS:
   INT sc_cmp (VEC_SCALAR x, const VEC_SCALAR y, const VECDATA_DESC *theVD);

   PARAMETERS:
   .  x - DOUBLE for each component of a vector data descriptor
   .  y - DOUBLE for each component of a vector data descriptor
   .  theVD - vector data descriptor

   DESCRIPTION:
   This function compares VEC_SCALARs.

   RETURN VALUE:
   INT
   .n    0 if VEC_SCALAR1 >  VEC_SCALAR2
   .n    1 if VEC_SCALAR1 <= VEC_SCALAR2
   D*/
/****************************************************************************/

INT sc_cmp (VEC_SCALAR x, const VEC_SCALAR y, const VECDATA_DESC *theVD)
{
  INT i;

  for (i=0; i<VD_NCOMP(theVD); i++)
    if (ABS(x[i])>=ABS(y[i]))
      return (0);

  return (1);
}

/****************************************************************************/
/*D
   sc_mul - x[i] = y[i] * z[i]

   SYNOPSIS:
   INT sc_mul (VEC_SCALAR x, const VEC_SCALAR y, VEC_SCALAR z, const VECDATA_DESC *theVD);

   PARAMETERS:
   .  x - DOUBLE for each component of a vector data descriptor
   .  y - DOUBLE for each component of a vector data descriptor
   .  z - DOUBLE for each component of a vector data descriptor
   .  theVD - vector data descriptor

   DESCRIPTION:
   This function calculates x[i] = y[i] * z[i] for every component
   of the 'VEC_SCALAR'.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    else if error occured.
   D*/
/****************************************************************************/

INT sc_mul (VEC_SCALAR x, const VEC_SCALAR y, const VEC_SCALAR z, const VECDATA_DESC *theVD)
{
  INT i;

  for (i=0; i<VD_NCOMP(theVD); i++)
    x[i] = y[i] * z[i];

  return (NUM_OK);
}

INT sc_mul_check (VEC_SCALAR x, const VEC_SCALAR y, const VEC_SCALAR z, const VECDATA_DESC *theVD)
{
  INT i;

  for (i=0; i<VD_NCOMP(theVD); i++)
  {
    x[i] = y[i] * z[i];
    if (x[i] == 0.0) x[i] = z[i];
  }
  return (NUM_OK);
}

/****************************************************************************/
/*D
   sc_read - Read VEC_SCALAR from input

   SYNOPSIS:
   INT sc_read (VEC_SCALAR x, const VECDATA_DESC *theVD,
   const char *name, INT argc, char **argv);

   PARAMETERS:
   .  x - DOUBLE for each component of a vector data descriptor
   .  theVD - vector data descriptor
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads VEC_SCALAR from input.
   It is used to read the arguments in 'npinit', e. g. the
   damping factors in the smoothers.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    else if error occured.
   D*/
/****************************************************************************/

#define OPTIONLEN                       32
#define OPTIONLENSTR            "31"
#define VALUELEN                        64
#define VALUELENSTR                     "63"

INT sc_read (VEC_SCALAR x, const VECDATA_DESC *theVD, const char *name, INT argc, char **argv)
{
  char option[OPTIONLEN],value[VALUELEN];
  INT i, n, found, type, err;
  const SHORT *offset;
  INT nDOUBLEs[NVECTYPES];
  DOUBLE theDOUBLEs[MAX_VEC_COMP][NVECTYPES];
  double lfValue;

  if (theVD != NULL)
    offset = VD_OFFSETPTR(theVD);
  if (strlen(name)>=OPTIONLEN-1) REP_ERR_RETURN (1);

  /* find input string */
  found = FALSE;
  for (i=0; i<argc; i++)
    if (sscanf(argv[i],expandfmt(CONCAT5("%",OPTIONLENSTR,"[a-zA-Z0-9_] %",VALUELENSTR,"[ -~]")),option,value)==2)
      if (strstr(option,name)!=NULL)
      {
        found = TRUE;
        break;
      }
  if (!found) REP_ERR_RETURN (2);

  /* read from value string */
  err = ReadVecTypeDOUBLEs(value,MAX_VEC_COMP,nDOUBLEs,theDOUBLEs);
  if (err!=NUM_OK)
    if (err==NUM_TYPE_MISSING)
    {
      /* iff no type is specified in the value string, scan one value for all */
      if (sscanf(value,"%lf",&lfValue)!=1)
        REP_ERR_RETURN (3);
      for (n=0; n<MAX_VEC_COMP; n++)
        x[n] = lfValue;
      return (NUM_OK);
    }
    else
      REP_ERR_RETURN (NUM_ERROR);

  /* fill x and check consistency with VECDATA_DESC */
  for (n=0, type=0; type<NVECTYPES; type++)
  {
    if (theVD!=NULL)
      if (n!=offset[type])
        REP_ERR_RETURN (4);
    for (i=0; i<nDOUBLEs[type]; i++)
      x[n++] = theDOUBLEs[i][type];
  }
  if (theVD!=NULL)
    if (n!=offset[type])
      REP_ERR_RETURN (4);

  return (NUM_OK);
}

/****************************************************************************/
/*D
   sc_disp - Display VEC_SCALAR

   SYNOPSIS:
   INT sc_disp (VEC_SCALAR x, const VECDATA_DESC *theVD, const char *name);

   PARAMETERS:
   .  x - DOUBLE for each component of a vector data descriptor
   .  theVD - vector data descriptor
   .  name - name of the argument

   DESCRIPTION:
   This function displays x on the shell.
   It is used to print the values of a VEC_SCALAR in 'npdisplay'.

   RETURN VALUE:
   INT
   .n    NUM_OK if ok
   .n    else if error occured.
   D*/
/****************************************************************************/

INT sc_disp (VEC_SCALAR x, const VECDATA_DESC *theVD, const char *name)
{
  char buffer[64];
  INT i, n, j, k;
  const SHORT *offset;

  UserWriteF(DISPLAY_NP_FORMAT_S,name);
  n = 0;
  if (theVD == NULL) {
    for (i=0; i<MAX_VEC_COMP; i++)
      if (i) UserWriteF("%s%-.4g",":",(double)x[n++]);
      else UserWriteF("%-.4g",(double)x[n++]);
    UserWrite("\n");
    return (NUM_OK);
  }

  offset = VD_OFFSETPTR(theVD);
  for (k=NVECTYPES; k>0; k--)
    if (offset[k]!=offset[k-1])
      break;

  for (i=0; i<k; i++)
  {
    if (i) UserWrite("|");
    switch (i)
    {
    case NODEVECTOR : UserWrite("nd "); break;
    case EDGEVECTOR : UserWrite("ed "); break;
    case ELEMVECTOR : UserWrite("el "); break;
    case SIDEVECTOR : UserWrite("sd "); break;
    }
    for (j=0; j<offset[i+1]-offset[i]; j++)
    {
      if (j) sprintf (buffer,"%s%-.4g",":",(double)x[n++]);
      else sprintf (buffer,"%-.4g",(double)x[n++]);
      UserWrite(buffer);
    }
  }
  UserWrite("\n");

  return (NUM_OK);
}
