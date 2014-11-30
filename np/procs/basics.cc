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

#include <config.h>
#include <cstring>

#include "ugdevices.h"

#include "numproc.h"
#include "np.h"
#include "general.h"
#include "scan.h"

#include "basics.h"

USING_UG_NAMESPACES

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

typedef struct
{
  NP_BASE base;

  DOUBLE a,b;
  VECDATA_DESC *f,*g,*d;

} NP_LC_VEC;

typedef struct
{
  NP_BASE base;

  VEC_SCALAR sp;
  VECDATA_DESC *x;
  VECDATA_DESC *y;

} NP_SCP_VEC;

typedef struct
{
  NP_BASE base;

  DOUBLE a;
  VECDATA_DESC *f;

} NP_SCALE_VEC;

typedef struct
{
  NP_BASE base;

  VECDATA_DESC *x;
  DOUBLE min,max;
  INT skip;

} NP_RANDOM_VEC;

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
   .  $s~<from~sym> - name of the source vector descriptor
   .  $d~<to~sym> - name of the destination vector descriptor

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
    scpv - numproc scalar product

   DESCRIPTION:
   This num proc performs a scalar product of two vectors

   'npinit <name>  $x <sym> $y <sym>'

   .  <name> - num proc name
   .  $x~<sym> - name of the first vector descriptor
   .  $y~<sym> - name of the second vector descriptor

   'npexecute <name>'

   EXAMPLE:
   .vb
   npcreate scpv $c copyv;
   npinit scpv $a a $b b;
   npexecute scpv;
   .ve
   D*/
/****************************************************************************/

static INT SCPV_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_SCP_VEC *np;

  np= (NP_SCP_VEC*)theNP;

  np->x = ReadArgvVecDesc(theNP->mg,"x",argc,argv);
  np->y = ReadArgvVecDesc(theNP->mg,"y",argc,argv);

  if ((np->x == NULL) || (np->y == NULL))
    return (NP_NOT_ACTIVE);

  return (NP_EXECUTABLE);
}

static INT SCPV_Display (NP_BASE *theNP)
{
  NP_SCP_VEC *np;

  np = (NP_SCP_VEC*)theNP;

  UserWrite("symbolic user data:\n");
  if (np->x != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(np->x));
  if (np->y != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"y",ENVITEM_NAME(np->y));
  sc_disp(np->sp,np->x,"scp");

  return (0);
}

static INT SCPV_Execute (NP_BASE *theNP, INT argc, char **argv)
{
  NP_SCP_VEC *np;

  np = (NP_SCP_VEC*)theNP;

  if ((np->x == NULL)  || (np->y == NULL)) return(1);

  if (ddotx(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ON_SURFACE,np->x,np->y,np->sp)) return (1);
  sc_disp(np->sp,np->x,"sp");

  return (0);
}

static INT SCPV_Construct (NP_BASE *theNP)
{
  theNP->Init = SCPV_Init;
  theNP->Display = SCPV_Display;
  theNP->Execute = SCPV_Execute;

  return(0);
}

/****************************************************************************/
/*D
   lcv - linear combination of two vectors

   DESCRIPTION:
   This num proc does a linear combination of two vectors

   'npinit <name>  $f <sym1> $g <sym2> $d <sym3> $a <value> $b <value>'

   .  <name> - num proc name
   .  $a~<value> - scaling factor of the first vector
   .  $f~<sym1> - name of the first source vector descriptor
   .  $b~<value> - scaling factor of the second vector
   .  $g~<sym2> - name of the second source vector descriptor
   .  $d~<sym3> - name of the destination vector descriptor

   'npexecute <name>'

   EXAMPLE:
   .vb
   npcreate lincomb $c lcv;
   npinit lincomb $a 1.0 $f sol1 $b -1.0 $g sol2 $d diff;
   npexecute lincomb;
   .ve
   D*/
/****************************************************************************/

static INT LCV_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_LC_VEC *np;

  np= (NP_LC_VEC*)theNP;

  np->f = ReadArgvVecDesc(theNP->mg,"f",argc,argv);
  np->g = ReadArgvVecDesc(theNP->mg,"g",argc,argv);
  np->d = ReadArgvVecDesc(theNP->mg,"d",argc,argv);
  if (np->d==NULL) np->d=np->f;
  if (ReadArgvDOUBLE("a",&np->a,argc,argv)) np->a = 1.0;
  if (ReadArgvDOUBLE("b",&np->b,argc,argv)) np->b = -1.0;

  if ((np->f == NULL) || (np->g == NULL)) return (NP_NOT_ACTIVE);
  return (NP_EXECUTABLE);
}

static INT LCV_Display (NP_BASE *theNP)
{
  NP_LC_VEC *np;

  np= (NP_LC_VEC*)theNP;

  UserWrite("symbolic user data:\n");
  if (np->f != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"f",ENVITEM_NAME(np->f));
  if (np->g != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"g",ENVITEM_NAME(np->g));
  if (np->d != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"d",ENVITEM_NAME(np->d));
  UserWriteF(DISPLAY_NP_FORMAT_SF,"a",(float)np->a);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"b",(float)np->b);

  return (0);
}

static INT LCV_Execute (NP_BASE *theNP, INT argc, char **argv)
{
  NP_LC_VEC *np;

  np= (NP_LC_VEC*)theNP;

  if ((np->f==NULL)  || (np->g==NULL) || (np->d==NULL) || np->f==np->g ) return(1);

  if (np->d!=np->f && np->d!=np->g)
  {
    if (dcopy(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,np->d,np->f)) return (1);
    if (dscal(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,np->d,np->a)) return (1);
    if (daxpy(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,np->d,np->b,np->g)) return (1);
  }
  if (np->d==np->f)
  {
    if (dscal(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,np->d,np->a)) return (1);
    if (daxpy(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,np->d,np->b,np->g)) return (1);
  }
  if (np->d==np->g)
  {
    if (dscal(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,np->d,np->b)) return (1);
    if (daxpy(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,np->d,np->a,np->f)) return (1);
  }

  return (0);
}

static INT LCV_Construct (NP_BASE *theNP)
{
  theNP->Init = LCV_Init;
  theNP->Display = LCV_Display;
  theNP->Execute = LCV_Execute;

  return(0);
}

/****************************************************************************/
/*D
   scalev - scale a vector

   DESCRIPTION:
   This num proc scales a vector by a factor

   'npinit <name>  $f <sym> $a <value>'

   .  <name> - num proc name
   .  $a~<value> - scaling factor
   .  $f~<sym> - name of the vector descriptor

   'npexecute <name>'

   EXAMPLE:
   .vb
   npcreate scale $c scalev;
   npinit scale $a -1.0 $f sol;
   npexecute scale;
   .ve
   D*/
/****************************************************************************/

static INT SCALEV_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_SCALE_VEC *np;

  np= (NP_SCALE_VEC*)theNP;

  np->f = ReadArgvVecDesc(theNP->mg,"f",argc,argv);
  if (ReadArgvDOUBLE("a",&np->a,argc,argv)) np->a = 1.0;

  if (np->f == NULL) return (NP_NOT_ACTIVE);
  return (NP_EXECUTABLE);
}

static INT SCALEV_Display (NP_BASE *theNP)
{
  NP_SCALE_VEC *np;

  np= (NP_SCALE_VEC*)theNP;

  UserWrite("symbolic user data:\n");
  if (np->f != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"f",ENVITEM_NAME(np->f));
  UserWriteF(DISPLAY_NP_FORMAT_SF,"a",(float)np->a);

  return (0);
}

static INT SCALEV_Execute (NP_BASE *theNP, INT argc, char **argv)
{
  NP_SCALE_VEC *np;

  np= (NP_SCALE_VEC*)theNP;

  if (np->f==NULL) return(1);
  if (dscal(NP_MG(theNP),0,CURRENTLEVEL(NP_MG(theNP)),ALL_VECTORS,np->f,np->a)) return (1);

  return (0);
}

static INT SCALEV_Construct (NP_BASE *theNP)
{
  theNP->Init = SCALEV_Init;
  theNP->Display = SCALEV_Display;
  theNP->Execute = SCALEV_Execute;

  return(0);
}

/****************************************************************************/
/*D
   rv - numproc random vector

   DESCRIPTION:
   The random vector numproc sets a vector to random values.

   'npinit <name> $x <vec sym> $min <min value> $max <max value> $skip [0|1]'

   .  <name> - num proc name
   .  $x~<vec~sym> - name of a vector symbol
   .  $min~<min~value> - minimal value, default 0
   .  $max~<max~value> - maximum value,default  1
   .  $skip~<skip> - 1 to skip SKIP-vectors,default 0

   'npexecute <name>'

   EXAMPLE:
   .vb
   npcreate rand $c rv;
   npinit rand $x x $min -1 $max 1 $skip 0;

   npexecute rand;
   .ve
   D*/
/****************************************************************************/

static INT RV_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_RANDOM_VEC *np;

  np = (NP_RANDOM_VEC*)theNP;
  np->x = ReadArgvVecDesc(theNP->mg,"x",argc,argv);
  if (np->x == NULL) return (NP_NOT_ACTIVE);
  if (ReadArgvDOUBLE("min",&np->min,argc,argv)) np->min = 0.0;
  if (ReadArgvDOUBLE("max",&np->max,argc,argv)) np->max = 1.0;
  if (ReadArgvINT("skip",&np->skip,argc,argv)) np->skip = 0;

  return (NP_EXECUTABLE);
}

static INT RV_Display (NP_BASE *theNP)
{
  NP_RANDOM_VEC *theRV;

  theRV   = (NP_RANDOM_VEC*)theNP;
  UserWrite("symbolic user data:\n");
  if (theRV->x != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(theRV->x));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"x","---");
  UserWriteF(DISPLAY_NP_FORMAT_SF,"min",(float)theRV->min);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"max",(float)theRV->max);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"skip",(float)theRV->skip);

  return (0);
}

static INT RV_Execute (NP_BASE *theNP, INT argc, char **argv)
{
  NP_RANDOM_VEC *theRV;
  GRID *g;
  INT i;

  theRV   = (NP_RANDOM_VEC*)theNP;
  if (theRV->x == NULL) return(1);

  for (i=0; i<=CURRENTLEVEL(NP_MG(theNP)); i++)
  {
    g=GRID_ON_LEVEL(NP_MG(theNP),i);
    if (l_dsetrandom2(g,theRV->x,EVERY_CLASS,theRV->min,theRV->max,theRV->skip)) return (1);
  }

  return (0);
}

static INT RV_Construct (NP_BASE *theNP)
{
  theNP->Init = RV_Init;
  theNP->Display = RV_Display;
  theNP->Execute = RV_Execute;

  return(0);
}

/****************************************************************************/
/** \brief Enrol basics

   This function creates the numprocs 'cv', 'cm' and 'eu'.
   It is called in InitNumerics.

   \return <ul>
   <li> 0 if ok </li>
   <li> 1 if error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX InitBasics ()
{
  if (CreateClass(BASE_CLASS_NAME ".cv",sizeof(NP_CLEAR_VEC),CV_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".cm",sizeof(NP_CLEAR_MAT),CM_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".eu",sizeof(NP_EUNORM_VEC),EU_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".copyv",sizeof(NP_COPY_VEC),COPYV_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".lcv",sizeof(NP_LC_VEC),LCV_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".scpv",sizeof(NP_SCP_VEC),SCPV_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".scalev",sizeof(NP_SCALE_VEC),SCALEV_Construct))
    return (__LINE__);
  if (CreateClass(BASE_CLASS_NAME ".rv",sizeof(NP_RANDOM_VEC),RV_Construct))
    return (__LINE__);

  return (0);
}
