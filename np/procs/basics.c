// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  basics.c                                                                                                      */
/*																			*/
/* Purpose:   basic numerical routines                                                                  */
/*																			*/
/* Author:	  Peter Bastian/Klaus Johannsen                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
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

#include "ugdevices.h"

#include "numproc.h"
#include "np.h"
#include "general.h"
#include "scan.h"

#include "basics.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define EU_STRDIR(p)                    ((p)->StructDir)

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct
{
  NP_BASE base;

  VECDATA_DESC *x;
  DOUBLE value;

} NP_CLEAR_VEC;

typedef struct
{
  NP_BASE base;

  MATDATA_DESC *A;
  DOUBLE value;

} NP_CLEAR_MAT;

typedef struct
{
  NP_BASE base;

  VECDATA_DESC *x;
  char StructDir[128];

} NP_EUNORM_VEC;

typedef struct
{
  NP_BASE base;

  VECDATA_DESC *s;
  VECDATA_DESC *d;

} NP_COPY_VEC;

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   cv - numproc clear vector

   DESCRIPTION:
   The clear vector numproc sets a vector to zero.

   'npinit <name> $x <vec sym>'

   .  <name> - num proc name
   .  $x~<vec~sym> - name of a vector symbol

   'npexecute <name>'

   EXAMPLE:
   .vb
   npcreate clearsol $c cv;
   npinit ckearsol $x x;

   npexecute clearsol;
   .ve
   D*/
/****************************************************************************/

static INT CV_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_CLEAR_VEC *np;

  np = (NP_CLEAR_VEC*)theNP;

  np->x = ReadArgvVecDesc(theNP->mg,"x",argc,argv);
  if (np->x == NULL)
    return (NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("value",&np->value,argc,argv))
    np->value = 0.0;

  return (NP_EXECUTABLE);
}

static INT CV_Display (NP_BASE *theNP)
{
  NP_CLEAR_VEC *theCV;

  theCV   = (NP_CLEAR_VEC*)theNP;
  UserWrite("symbolic user data:\n");
  if (theCV->x != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(theCV->x));
  UserWriteF(DISPLAY_NP_FORMAT_SF,"value",(float)theCV->value);

  return (0);
}

static INT CV_Execute (NP_BASE *theNP, INT argc, char **argv)
{
  NP_CLEAR_VEC *theCV;
  INT cl;

  theCV   = (NP_CLEAR_VEC*)theNP;
  if (theCV->x == NULL) return(1);

  cl = CURRENTLEVEL(NP_MG(theNP));
  if (dset(NP_MG(theNP),0,cl,ALL_VECTORS,theCV->x,theCV->value))
    return (1);

  return (0);
}

static INT CV_Construct (NP_BASE *theNP)
{
  theNP->Init = CV_Init;
  theNP->Display = CV_Display;
  theNP->Execute = CV_Execute;

  return(0);
}

/****************************************************************************/
/*D
   cm - numproc clear matrix

   DESCRIPTION:
   The clear matrix numproc sets a matrix to zero.

   'npinit <name> $A <mat sym>'

   .  <name> - num proc name
   .  $A~<mat~sym> - matrix descriptor

   'npexecute <name>'

   EXAMPLE:
   .vb
   npcreate clearA $c cm;
   npinit clearA $A MAT;

   npexecute clearA;
   .ve
   D*/
/****************************************************************************/

static INT CM_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_CLEAR_MAT *np;

  np = (NP_CLEAR_MAT*)theNP;

  np->A = ReadArgvMatDesc(theNP->mg,"A",argc,argv);
  if (np->A == NULL)
    return (NP_NOT_ACTIVE);

  if (ReadArgvDOUBLE("value",&np->value,argc,argv))
    np->value = 0.0;

  return (NP_EXECUTABLE);
}

static INT CM_Display (NP_BASE *theNP)
{
  NP_CLEAR_MAT *theCM;

  theCM   = (NP_CLEAR_MAT*)theNP;

  UserWrite("symbolic user data:\n");
  if (theCM->A != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(theCM->A));
  UserWriteF(DISPLAY_NP_FORMAT_SF,"value",(float)theCM->value);

  return (0);
}

static INT CM_Execute (NP_BASE *theNP, INT argc, char **argv)
{
  NP_CLEAR_MAT *theCM;

  theCM   = (NP_CLEAR_MAT*)theNP;
  if (theCM->A == NULL) return(1);

  if (dmatset(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,theCM->A,theCM->value))
    return (1);

  return (0);
}

static INT CM_Construct (NP_BASE *theNP)
{
  theNP->Init = CM_Init;
  theNP->Display = CM_Display;
  theNP->Execute = CM_Execute;

  return(0);
}

/****************************************************************************/
/*D
   eu - numproc to compute euklidian norm of a vector

   DESCRIPTION:
   The numproc computes the euklidian norm of a vector.

   'npinit <name> $x <vec sym>'

   .  <name> - num proc name
   .  $x~<vec~sym> - name of a vector descriptor

   'npexecute <name>'

   EXAMPLE:
   .vb
   npcreate euklid $c eu;
   npinit euklid $x cor;

   npexecute euklid;
   .ve
   D*/
/****************************************************************************/

static INT EU_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_EUNORM_VEC *theEU;
  INT i;

  theEU   = (NP_EUNORM_VEC*)theNP;
  EU_STRDIR(theEU)[0] = '\0';
  for (i=1; i<argc; i++)
    if (argv[i][0]=='s')
    {
      if (sscanf(argv[i],"s %s",EU_STRDIR(theEU))!=1)
        EU_STRDIR(theEU)[0] = '\0';
      break;
    }

  theEU->x = ReadArgvVecDesc(theNP->mg,"x",argc,argv);
  if (theEU->x == NULL)
    return (NP_NOT_ACTIVE);

  return (NP_EXECUTABLE);
}

static INT EU_Display (NP_BASE *theNP)
{
  NP_EUNORM_VEC *theEU;

  theEU = (NP_EUNORM_VEC*)theNP;

  UserWrite("symbolic user data:\n");
  if (theEU->x != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(theEU->x));
  UserWrite("configuration parameters:\n");
  if (EU_STRDIR(theEU)[0]!='\0')
    UserWriteF(DISPLAY_NP_FORMAT_SS,"structdir",EU_STRDIR(theEU));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"structdir","---");

  return (0);
}

static INT EU_Execute (NP_BASE *theNP, INT argc, char **argv)
{
  NP_EUNORM_VEC *theEU;
  VEC_SCALAR eunorm;
  INT cl;

  theEU   = (NP_EUNORM_VEC*)theNP;

  if (theEU->x == NULL) return(1);

  cl = CURRENTLEVEL(NP_MG(theNP));
  if (dnrm2x(NP_MG(theNP),cl,cl,ALL_VECTORS,theEU->x,eunorm))
    return (1);

  if (WriteVEC_SCALAR(theEU->x,eunorm,EU_STRDIR(theEU)))
    return (1);

  return (0);
}

static INT EU_Construct (NP_BASE *theNP)
{
  theNP->Init = EU_Init;
  theNP->Display = EU_Display;
  theNP->Execute = EU_Execute;

  return(0);
}

/****************************************************************************/
/*D
   copyv - numproc copy vector

   DESCRIPTION:
   This num proc copys a vector

   'npinit <name>  $f <from sym> $t <to sym>'

   .  <name> - num proc name
   .  $f~<from~sym> - name of the source vector descriptor
   .  $t~<to~sym> - name of the destination vector descriptor

   'npexecute <name>'

   EXAMPLE:
   .vb
   npcreate copya2b $c copyv;
   npinit copya2b $s a $d b;
   npexecute copya2b;
   .ve
   D*/
/****************************************************************************/

static INT COPYV_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_COPY_VEC *np;

  np= (NP_COPY_VEC*)theNP;

  np->s = ReadArgvVecDesc(theNP->mg,"s",argc,argv);
  np->d = ReadArgvVecDesc(theNP->mg,"d",argc,argv);

  if ((np->s == NULL) || (np->d == NULL))
    return (NP_NOT_ACTIVE);

  return (NP_EXECUTABLE);
}

static INT COPYV_Display (NP_BASE *theNP)
{
  NP_COPY_VEC *theCopyV;

  theCopyV = (NP_COPY_VEC*)theNP;

  UserWrite("symbolic user data:\n");
  if (theCopyV->s != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"s",ENVITEM_NAME(theCopyV->s));
  if (theCopyV->d != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"d",ENVITEM_NAME(theCopyV->d));

  return (0);
}

static INT COPYV_Execute (NP_BASE *theNP, INT argc, char **argv)
{
  NP_COPY_VEC *np;

  np = (NP_COPY_VEC*)theNP;

  if ((np->d == NULL)  || (np->s == NULL)) return(1);

  if (dcopy(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,np->d,np->s))
    return (1);

  return (0);
}

static INT COPYV_Construct (NP_BASE *theNP)
{
  theNP->Init = COPYV_Init;
  theNP->Display = COPYV_Display;
  theNP->Execute = COPYV_Execute;

  return(0);
}

/****************************************************************************/
/*D
   InitBasics - Enrol basics

   SYNOPSIS:
   INT InitBasics (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function creates the numprocs 'cv', 'cm' and 'eu'.
   It is called in InitNumerics.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT InitBasics (void)
{
  if (CreateClass(BASE_CLASS_NAME ".cv",sizeof(NP_CLEAR_VEC),CV_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".cm",sizeof(NP_CLEAR_MAT),CM_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".eu",sizeof(NP_EUNORM_VEC),EU_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".copyv",sizeof(NP_COPY_VEC),
                  COPYV_Construct))
    return (__LINE__);

  return (0);
}
