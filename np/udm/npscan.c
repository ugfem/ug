// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  npscan.c	                                                                                                */
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
/*			  np part of former np/udm/scan.c, 15.5.97						*/
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
#include <ctype.h>
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
#define BLANKS                  " \t\n"

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

  return(1);
}

/****************************************************************************/
/*D
   ReadArgvVecDescX - Read command strings

   SYNOPSIS:
   VECDATA_DESC *ReadArgvVecDescX (MULTIGRID *theMG, const char *name,
                                                           INT argc, char **argv, INT CreateIfNonExistent);

   PARAMETERS:
   .  theMG - pointer to a multigrid
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector
   .  CreateIfNonExistent - create vector if not yet existing

   DESCRIPTION:
   This function reads a symbol name from the command strings and returns
   a pointer to the corresponding vector descriptor.

   CAUTION: If no template is specified the first vector template is used.

   This call locks the vector descriptor for dynamic allocation.

   SYNTAX:
   '$<name> <vec desc name>[/<template name>]'

   RETURN VALUE:
   VECDATA_DESC *
   .n    pointer to vector descriptor
   .n    NULL if error occurs
   D*/
/****************************************************************************/

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
   (This function is a macro calling 'ReadArgvVecDescX')

   CAUTION: If no template is specified the first vector template is used.

   This call locks the vector descriptor for dynamic allocation.

   SYNTAX:
   '$<name> <vec desc name>[/<template name>]'

   RETURN VALUE:
   VECDATA_DESC *
   .n    pointer to vector descriptor
   .n    NULL if error occurs
   D*/
/****************************************************************************/

VECDATA_DESC *ReadArgvVecDescX (MULTIGRID *theMG, const char *name,
                                INT argc, char **argv, INT CreateIfNonExistent)
{
  VECDATA_DESC *vd;
  char value[VALUELEN],vdname[NAMESIZE],tname[NAMESIZE];
  INT res;

  ASSERT(theMG!=NULL);

  if (ReadArgvChar(name,value,argc,argv))
    return(NULL);

  res = sscanf(value,expandfmt(CONCAT5("%",NAMELENSTR,"[a-zA-Z0-9_] / %",NAMELENSTR,"[a-zA-Z0-9_]")),vdname,tname);
  vd = GetVecDataDescByName(theMG,vdname);
  if (vd == NULL && CreateIfNonExistent)
  {
    if (res==2)
      vd = CreateVecDescOfTemplate (theMG,vdname,tname);
    else
      /* taking default template */
      vd = CreateVecDescOfTemplate (theMG,vdname,NULL);
  }
  if (vd == NULL) return(NULL);

  if (LockVD(vd)) REP_ERR_RETURN (NULL);

  return(vd);
}

/****************************************************************************/
/*D
   ReadArgvVecTemplate - read vec template from command string

   SYNOPSIS:
   VEC_TEMPLATE *ReadArgvVecTemplate (const FORMAT *fmt, const char *name,
                                                           INT argc, char **argv)

   PARAMETERS:
   .  fmt - pointer to a format
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector

   DESCRIPTION:
   This function reads a vector template name from the command strings and returns
   a pointer to the corresponding vector template.

   SYNTAX:
   '$<name> <vec template name>'

   RETURN VALUE:
   VECDATA_DESC *
   .n    pointer to vector descriptor
   .n    NULL if error occurs
   D*/
/****************************************************************************/

VEC_TEMPLATE *ReadArgvVecTemplate (const FORMAT *fmt, const char *name,
                                   INT argc, char **argv)
{
  char value[VALUELEN],vtname[NAMESIZE];

  if (ReadArgvChar(name,value,argc,argv))
    return(NULL);

  if (sscanf(value,expandfmt(CONCAT3("%",NAMELENSTR,"[a-zA-Z0-9_]")),vtname)!=1)
    return(NULL);

  return(GetVectorTemplate(fmt,vtname));
}

/****************************************************************************/
/*D
   ReadArgvVecTemplateSub - read vec template sub descriptor from command string

   SYNOPSIS:
   VEC_TEMPLATE *ReadArgvVecTemplateSub (const FORMAT *fmt, const char *name,
                                                                           INT argc, char **argv, INT *sub)

   PARAMETERS:
   .  fmt - pointer to a format
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector
   .  sub - index of sub descriptor in template

   DESCRIPTION:
   This function reads a vector template name and the name of one of its sub descriptors
   from the command strings and returns a pointer to the corresponding vector template plus
   the index of the sub descriptor in the template.
   If the sub descriptor is omitted the template itself is meant and sub will contain FULL_TPLT.

   SYNTAX:
   '$<name> <vec template name> [<sub descriptor name>]'

   RETURN VALUE:
   VECDATA_DESC *
   .n    pointer to vector descriptor
   .n    NULL if error occurs
   D*/
/****************************************************************************/

VEC_TEMPLATE *ReadArgvVecTemplateSub (const FORMAT *fmt, const char *name,
                                      INT argc, char **argv, INT *sub)
{
  VEC_TEMPLATE *vt;
  INT i,res;
  char value[VALUELEN],vtname[NAMESIZE],subname[NAMESIZE];

  if (ReadArgvChar(name,value,argc,argv))
    return(NULL);

  res = sscanf(value,expandfmt(CONCAT5("%",NAMELENSTR,"[a-zA-Z0-9_] %",NAMELENSTR,"[a-zA-Z0-9_]")),vtname,subname);
  if (res<1)
    return (NULL);

  vt = GetVectorTemplate(fmt,vtname);
  if (vt==NULL)
    return(NULL);

  if (res==1)
  {
    *sub = FULL_TPLT;
    return (vt);
  }

  for (i=0; i<VT_NSUB(vt); i++)
    if (strcmp(SUBV_NAME(VT_SUB(vt,i)),subname)==0)
      break;
  if (i>=VT_NSUB(vt))
    return(NULL);

  *sub = i;

  return (vt);
}

/****************************************************************************/
/*D
   ReadArgvMatTemplateSub - read vec template sub descriptor from command string

   SYNOPSIS:
   MAT_TEMPLATE *ReadArgvMatTemplateSub (const FORMAT *fmt, const char *name,
                                                                           INT argc, char **argv, INT *sub)

   PARAMETERS:
   .  fmt - pointer to a format
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector
   .  sub - index of sub descriptor in template

   DESCRIPTION:
   This function reads a vector template name and the name of one of its sub descriptors
   from the command strings and returns a pointer to the corresponding vector template plus
   the index of the sub descriptor in the template.
   If the sub descriptor is omitted the template itself is meant and sub will contain FULL_TPLT.

   SYNTAX:
   '$<name> <vec template name> <sub descriptor name>'

   RETURN VALUE:
   VECDATA_DESC *
   .n    pointer to vector descriptor
   .n    NULL if error occurs
   D*/
/****************************************************************************/

MAT_TEMPLATE *ReadArgvMatTemplateSub (const FORMAT *fmt, const char *name,
                                      INT argc, char **argv, INT *sub)
{
  MAT_TEMPLATE *mt;
  INT i,res;
  char value[VALUELEN],mtname[NAMESIZE],subname[NAMESIZE];

  if (ReadArgvChar(name,value,argc,argv))
    return (NULL);

  res = sscanf(value,expandfmt(CONCAT5("%",NAMELENSTR,"[a-zA-Z0-9_] %",NAMELENSTR,"[a-zA-Z0-9_]")),mtname,subname);
  if (res<1)
    return (NULL);

  mt = GetMatrixTemplate(fmt,mtname);
  if (mt==NULL)
    return(NULL);

  if (res==1)
  {
    *sub = FULL_TPLT;
    return (mt);
  }

  for (i=0; i<MT_NSUB(mt); i++)
    if (strcmp(SUBM_NAME(MT_SUB(mt,i)),subname)==0)
      break;
  if (i>=MT_NSUB(mt))
    return(NULL);

  *sub = i;

  return (mt);
}

/****************************************************************************/
/*D
   ReadArgvMatDescX - Read command strings

   SYNOPSIS:
   MATDATA_DESC *ReadArgvMatDescX (MULTIGRID *theMG, const char *name,
                                                           INT argc, char **argv, INT CreateIfNonExistent);

   PARAMETERS:
   .  theMG - pointer to a multigrid
   .  name - name of the argument
   .  argc - argument counter
   .  argv - argument vector
   .  CreateIfNonExistent - create vector if not yet existing

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
   (This function is a macro calling 'ReadArgvMatDescX')

   RETURN VALUE:
   MATDATA_DESC *
   .n    pointer to matrix descriptor
   .n    NULL if error occurs
   D*/
/****************************************************************************/

MATDATA_DESC *ReadArgvMatDescX (MULTIGRID *theMG, const char *name,
                                INT argc, char **argv, INT CreateIfNonExistent)
{
  MATDATA_DESC *md;
  char value[VALUELEN],mdname[NAMESIZE],tname[NAMESIZE];
  INT res;

  ASSERT(theMG!=NULL);

  if (ReadArgvChar(name,value,argc,argv))
    return (NULL);

  res = sscanf(value,expandfmt(CONCAT5("%",NAMELENSTR,
                                       "[a-zA-Z0-9_] / %",
                                       NAMELENSTR,"[a-zA-Z0-9_]")),
               mdname,tname);
  md = GetMatDataDescByName(theMG,mdname);
  if (md==NULL && CreateIfNonExistent)
  {
    if (res==2)
      md = CreateMatDescOfTemplate (theMG,mdname,tname);
    else
      /* taking default template */
      md = CreateMatDescOfTemplate (theMG,mdname,NULL);
  }
  if (md == NULL) return (NULL);

  if (LockMD(md)) REP_ERR_RETURN(NULL);

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

  ASSERT(theMG!=NULL);

  if (ReadArgvChar(name,value,argc,argv))
    return (NULL);

  return(GetNumProcByName(theMG,value,class));
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
/*
   ReadVecTypeINTs - Read a number of INTs from the input string

   SYNOPSIS:
   INT ReadVecTypeINTs (const FORMAT *fmt, char *str, INT n, INT nINT[MAXVECTORS], INT theINTs[][MAXVECTORS]);

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

INT ReadVecTypeINTs (const FORMAT *fmt, char *str, INT n, INT nINT[MAXVECTORS], INT theINTs[][MAXVECTORS])
{
  char *s,*tok,*typetok[MAXVECTORS];
  INT type;
  int iValue;

  for (type=0; type<MAXVECTORS; type++) {nINT[type] = 0; typetok[type] = NULL;}

  for (tok=strtok(str,TYPESEP); tok!=NULL; tok=strtok(NULL,TYPESEP))
  {
    /* find first character indicating abstract vtype */
    for (s=tok; *s!='\0'; s++)
      if (strchr(BLANKS,*s)==NULL)
        break;
    if (!isalpha(*s) || ((type = FMT_N2T(fmt,*s))==NOVTYPE))
    {
      PrintErrorMessageF('E',"ReadVecTypeINTs","could not read type specifier or invalid type (in '%s')\n",str);
      REP_ERR_RETURN (1);
    }
    typetok[type] = ++s;
    if (isalpha(*s))
    {
      PrintErrorMessage('E',"ReadVecTypeINTs","two chars for vtype specification is not supported anymore\n"
                        "please read the CHANGES from ug-3.7 to ug-3.8");
      REP_ERR_RETURN (2);
    }
  }

  for (type=0; type<MAXVECTORS; type++)
    if (typetok[type]!=NULL)
      for (tok=strtok(typetok[type],COMPSEP); tok!=NULL; tok=strtok(NULL,COMPSEP))
      {
        if (nINT[type]>=n)
        {
          PrintErrorMessageF('E',"ReadVecTypeINTs","max number of INTs exceeded (in '%s')\n",str);
          REP_ERR_RETURN (3);
        }

        if (sscanf(tok,"%d",&iValue)!=1)
        {
          PrintErrorMessageF('E',"ReadVecTypeINTs","could not scan INT (in '%s')\n",str);
          return (3);
        }
        else
          theINTs[nINT[type]++][type] = (INT) iValue;
      }

  return (NUM_OK);
}

/****************************************************************************/
/*
   ReadVecTypeDOUBLEs - Read a number of DOUBLEs from the input string

   SYNOPSIS:
   INT ReadVecTypeDOUBLEs (const FORMAT *fmt, char *str, INT n,
                                        INT nDOUBLE[MAXVECTORS], DOUBLE theDOUBLEs[][MAXVECTORS]);

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

INT ReadVecTypeDOUBLEs (const FORMAT *fmt, char *str, INT n, INT nDOUBLE[MAXVECTORS], DOUBLE theDOUBLEs[][MAXVECTORS])
{
  char *s,*tok,*typetok[MAXVECTORS],*notypetok;
  INT type,found;
  double lfValue;

  for (type=0; type<MAXVECTORS; type++) {nDOUBLE[type] = 0; typetok[type] = NULL;}
  notypetok = NULL;

  for (tok=strtok(str,TYPESEP); tok!=NULL; tok=strtok(NULL,TYPESEP))
  {
    /* find first character indicating abstract vtype */
    for (s=tok; *s!='\0'; s++)
      if (strchr(BLANKS,*s)==NULL)
        break;
    if (!isalpha(*s) || ((type = FMT_N2T(fmt,*s))==NOVTYPE))
      notypetok = tok;
    else
    {
      typetok[type] = ++s;
      if (isalpha(*s))
      {
        PrintErrorMessage('E',"ReadVecTypeDOUBLEs","two chars for vtype specification is not supported anymore\n"
                          "please read the CHANGES from ug-3.7 to ug-3.8");
        REP_ERR_RETURN (2);
      }
    }
  }

  found = 0;
  for (type=0; type<MAXVECTORS; type++)
    if (typetok[type]!=NULL)
      for (tok=strtok(typetok[type],COMPSEP); tok!=NULL; tok=strtok(NULL,COMPSEP))
      {
        found++;
        if (nDOUBLE[type]>=n)
        {
          PrintErrorMessageF('E',"ReadVecTypeDOUBLEs","max number of DOUBLEs exceeded (in '%s')\n",str);
          REP_ERR_RETURN (3);
        }

        if (sscanf(tok,"%lf",&lfValue)!=1)
        {
          PrintErrorMessageF('E',"ReadVecTypeDOUBLEs","could not scan DOUBLE (in '%s')\n",str);
          return (3);
        }
        else
          theDOUBLEs[nDOUBLE[type]++][type] = (DOUBLE) lfValue;
      }

  if (notypetok!=NULL)
  {
    if (found)
    {
      PrintErrorMessageF('E',"ReadVecTypeDOUBLEs","type specifier missing (in '%s')\n",str);
      REP_ERR_RETURN (NUM_ERROR);
    }

    /* there is only one token without type label */

    /* is there only one value? */
    found = 0;
    for (tok=strtok(notypetok,COMPSEP); tok!=NULL; tok=strtok(NULL,COMPSEP))
      found++;
    if (found!=1)
    {
      PrintErrorMessageF('E',"ReadVecTypeDOUBLEs","type specifier missing but several values given (in '%s')\n",str);
      REP_ERR_RETURN (NUM_ERROR);
    }
    else
      return (NUM_TYPE_MISSING);
  }

  return (NUM_OK);
}

/****************************************************************************/
/*
   ReadVecTypeOrder - Read a number of INTs from the input string

   SYNOPSIS:
   INT ReadVecTypeOrder (const FORMAT *fmt, char *str, INT n,
                                INT MaxPerType, INT *nOrder, INT theOrder[]);

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

INT ReadVecTypeOrder (const FORMAT *fmt, char *str, INT n, INT MaxPerType, INT *nOrder, INT theOrder[])
{
  char *token,c;
  INT ni,tp;
  int iValue;

  ni = 0;
  for (token=strtok(str,COMPSEP); token!=NULL; token=strtok(NULL,COMPSEP))
  {
    if (ni>=n)
    {
      PrintErrorMessageF('E',"ReadVecTypeOrder","max number of values exceeded (in '%s')\n",str);
      REP_ERR_RETURN (3);
    }

    if              ((sscanf(token,"%c%d",&c,&iValue)==2) && (iValue<MaxPerType))
    {
      if ((tp=FMT_N2T(fmt,c))==NOVTYPE)
      {
        PrintErrorMessageF('E',"ReadVecTypeOrder","invalid type specified (in '%s')\n",str);
        REP_ERR_RETURN (2);
      }
      else
        theOrder[ni++] = tp*MaxPerType + (INT) iValue;
    }
    else
    {
      PrintErrorMessage('E',"ReadVecTypeOrder","two chars for vtype specification is not supported anymore\n"
                        "please read the CHANGES from ug-3.7 to ug-3.8");
      REP_ERR_RETURN (3);
    }
  }

  *nOrder = ni;

  return (NUM_OK);
}

/****************************************************************************/
/*
   ReadVecTypeNUMPROCs - Read a number of NUMPROCs from the input string

   SYNOPSIS:
   INT ReadVecTypeNUMPROCs (const MULTIGRID *theMG, char *str, char *class_name, INT n, INT nNUMPROC[MAXVECTORS], NP_BASE *theNUMPROCs[][MAXVECTORS]);

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

INT ReadVecTypeNUMPROCs (const MULTIGRID *theMG, char *str, char *class_name, INT n, INT nNUMPROC[MAXVECTORS], NP_BASE *theNUMPROCs[][MAXVECTORS])
{
  FORMAT *fmt;
  char *s,*tok,*typetok[MAXVECTORS];
  INT type;

  ASSERT(theMG!=NULL);

  for (type=0; type<MAXVECTORS; type++) {nNUMPROC[type] = 0; typetok[type] = NULL;}

  fmt = MGFORMAT(theMG);
  for (tok=strtok(str,TYPESEP); tok!=NULL; tok=strtok(NULL,TYPESEP))
  {
    /* find first character indicating abstract vtype */
    for (s=tok; *s!='\0'; s++)
      if (strchr(BLANKS,*s)==NULL)
        break;
    if (!isalpha(*s) || ((type = FMT_N2T(fmt,*s))==NOVTYPE))
    {
      PrintErrorMessageF('E',"ReadVecTypeNUMPROCs","could not read type specifier or invalid type (in '%s')\n",str);
      REP_ERR_RETURN (1);
    }
    typetok[type] = ++s;
    if (isalpha(*s))
    {
      PrintErrorMessage('E',"ReadVecTypeNUMPROCs","two chars for vtype specification is not supported anymore\n"
                        "please read the CHANGES from ug-3.7 to ug-3.8");
      REP_ERR_RETURN (2);
    }
  }

  for (type=0; type<MAXVECTORS; type++)
    if (typetok[type]!=NULL)
      for (tok=strtok(typetok[type],COMPSEP); tok!=NULL; tok=strtok(NULL,COMPSEP))
      {
        if (nNUMPROC[type]>=n)
        {
          PrintErrorMessageF('E',"ReadVecTypeNUMPROCs","max number of NUMPROCs exceeded (in '%s')\n",str);
          REP_ERR_RETURN (3);
        }

        if ((theNUMPROCs[nNUMPROC[type]++][type]=GetNumProcByName(theMG,tok,class_name))==NULL)
        {
          PrintErrorMessageF('E',"ReadVecTypeNUMPROCs","NUMPROC '%s' not found (in '%s')\n",tok,str);
          REP_ERR_RETURN (3);
        }
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

  if (VD_NID(theVD)!=NO_IDENT)
  {
    DOUBLE sx,sy;
    INT j;

    for (i=0; i<VD_NCOMP(theVD); i++)
      if (VD_IDENT(theVD,i)==i)
      {
        sx = sy = 0;
        for (j=0; j<VD_NCOMP(theVD); j++)
          if (VD_IDENT(theVD,j)==i)
          {
            sx += x[j]*x[j];
            sy += y[j]*y[j];
          }
        if (sx>=sy)
          return (0);
      }
    return (1);
  }

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
   INT sc_read (VEC_SCALAR x, const FORMAT *fmt, const VECDATA_DESC *theVD, const char *name, INT argc, char **argv);

   PARAMETERS:
   .  x - DOUBLE for each component of a vector data descriptor
   .  theVD - vector data descriptor (may be NULL)
   .  fmt - corresponding format (theVD my be NULL)
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

INT sc_read (VEC_SCALAR x, const FORMAT *fmt, const VECDATA_DESC *theVD, const char *name, INT argc, char **argv)
{
  char option[OPTIONLEN],value[VALUELEN];
  INT i, n, found, type, err;
  const SHORT *offset;
  INT nDOUBLEs[NVECTYPES];
  DOUBLE theDOUBLEs[MAX_VEC_COMP][NVECTYPES];
  double lfValue;

  if (theVD != NULL)
  {
    offset = VD_OFFSETPTR(theVD);
    if (fmt!=MGFORMAT(VD_MG(theVD)))
      REP_ERR_RETURN (1);
  }
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
  if (!found) return(2);

  /* read from value string */
  err = ReadVecTypeDOUBLEs(fmt,value,MAX_VEC_COMP,nDOUBLEs,theDOUBLEs);
  if (err!=NUM_OK)
    if (err==NUM_TYPE_MISSING)
    {
      /* iff no type is specified in the value string, scan one value for all */
      if (sscanf(value,"%lf",&lfValue)!=1)
      {
        PrintErrorMessageF('E',"sc_read","could not scan single value (in '%s')\n",value);
        REP_ERR_RETURN (3);
      }
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
      {
        PrintErrorMessageF('E',"sc_read","number of values per type does not coincide with vd (in '%s')\n",value);
        REP_ERR_RETURN (4);
      }
    for (i=0; i<nDOUBLEs[type]; i++)
      x[n++] = theDOUBLEs[i][type];
  }
  if (theVD!=NULL)
    if (n!=offset[type])
    {
      PrintErrorMessageF('E',"sc_read","total number of values does not coincide with vd (in '%s')\n",value);
      REP_ERR_RETURN (4);
    }

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
  const FORMAT *fmt;
  INT i, n, j, k;
  const SHORT *offset;
  char c;

  UserWriteF(DISPLAY_NP_FORMAT_S,name);
  n = 0;
  if (theVD == NULL) {
    for (i=0; i<MAX_VEC_COMP; i++)
      if (i) UserWriteF("%s%-.4g",":",(double)x[n++]);
      else UserWriteF("%-.4g",(double)x[n++]);
    UserWrite("\n");
    return (NUM_OK);
  }

  fmt = MGFORMAT(VD_MG(theVD));
  offset = VD_OFFSETPTR(theVD);
  for (k=NVECTYPES; k>0; k--)
    if (offset[k]!=offset[k-1])
      break;

  for (i=0; i<k; i++)
  {
    if (i) UserWrite(TYPESEP);
    c = FMT_T2N(fmt,i);
    UserWriteF("%c  ",c);
    for (j=0; j<offset[i+1]-offset[i]; j++)
    {
      if (j) UserWriteF("%s%-.4g",":",(double)x[n++]);
      else UserWriteF("%-.4g",(double)x[n++]);
    }
  }
  UserWrite("\n");

  return (NUM_OK);
}
