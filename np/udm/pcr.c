// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  pcr.c	                                                                                                        */
/*																			*/
/* Purpose:   print convergence rates                                                                           */
/*																			*/
/* Author:	  Henrik Rentz-Reichert/Klaus Johannsen                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   25.03.95 begin, ug version 3.0								*/
/*			  09.12.95 transition to new descriptor formats (HRR)			*/
/*			  December 3, 1996 new np subsystem                                                     */
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "compiler.h"
#include "ugdevices.h"
#include "misc.h"
#include "gm.h"
#include "algebra.h"
#include "evm.h"
#include "ugstruct.h"
#include "debug.h"
#include "general.h"
#include "np.h"

#include "pcr.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

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

/* general purpose text buffer */
static char buffer[256];

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* Function:  GetStrINTinRange												*/
/*																			*/
/* Purpose:   get INT value from str and check wether it is lying in the	*/
/*			  specified interval											*/
/*																			*/
/* Input:	  s.a															*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT GetStrINTinRange (const char *str, INT min, INT max, INT *value)
{
  int iValue;

  if (sscanf(str,"%d",&iValue)!=1)
  {
    PrintErrorMessageF('E',"GetStrINTinRange","could not scan INT value from string '%s'",str);
    return(2);
  }
  if (iValue<min)
  {
    PrintErrorMessageF('E',"GetStrINTinRange","value (%d) < min (%g)",iValue,min);
    return(3);
  }
  if (iValue>max)
  {
    PrintErrorMessageF('E',"GetStrINTinRange","value (%d) > max (%g)",iValue,max);
    return(4);
  }
  *value = (INT) iValue;

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  GetStrDOUBLEinRange											*/
/*																			*/
/* Purpose:   get DOUBLE value from str and check wether it is lying in the	*/
/*			  specified interval											*/
/*																			*/
/* Input:	  s.a															*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT GetStrDOUBLEinRange (const char *str, DOUBLE min, DOUBLE max, DOUBLE *value)
{
  float fValue;

  if (sscanf(str,"%f",&fValue)!=1)
  {
    PrintErrorMessageF('E',"GetStrDOUBLEinRange","could not scan DOUBLE value from string '%s'",str);
    return(2);
  }
  if (fValue<min)
  {
    PrintErrorMessageF('E',"GetStrDOUBLEinRange","value (%d) < min (%g)",fValue,min);
    return(3);
  }
  if (fValue>max)
  {
    PrintErrorMessageF('E',"GetStrDOUBLEinRange","value (%d) > max (%g)",fValue,max);
    return(4);
  }
  *value = (DOUBLE) fValue;

  return (0);
}

/****************************************************************************/
/*D
   WriteVEC_SCALAR - Write VEC_SCALAR on sreen and to stringvariables

   SYNOPSIS:
   INT WriteVEC_SCALAR (const VECDATA_DESC *theVDT, const VEC_SCALAR Scalar, const char *structdir);

   PARAMETERS:
   .  theVDT - type vector descriptor
   .  Scalar - DOUBLE for each component of a type vector descriptor
   .  structdir - ugstruct for stringvariables

   DESCRIPTION:
   This function writes VEC_SCALAR on the shell and to stringvariables.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT WriteVEC_SCALAR (const VECDATA_DESC *theVDT, const VEC_SCALAR Scalar, const char *structdir)
{
  INT i;
  char name[2];

  for (i=0; i<VD_NCOMP(theVDT); i++)
    UserWriteF("%c: %-12.7e\n",VM_COMP_NAME(theVDT,i),Scalar[i]);

  if (structdir[0]!='\0')
  {
    if (ChangeStructDir(structdir)==NULL) REP_ERR_RETURN (1);
    for (i=0; i<VD_NCOMP(theVDT); i++)
    {
      sprintf(name,"%c",VM_COMP_NAME(theVDT,i));
      if (SetStringValue(name,Scalar[i])) REP_ERR_RETURN (1);
    }
    if (ChangeStructDir(":")==NULL) REP_ERR_RETURN (1);
  }
  return (0);
}

/****************************************************************************/
/*
   PrintConvergenceRate - Print convergence rate

   SYNOPSIS:
   INT PreparePCR (VECDATA_DESC *theVDT, INT DispMode, const char *text,
   INT *ID);

   PARAMETERS:
   .  theVDT - vector data descriptor
   .  DispMode - display modus
   .  text - text printed on the shell
   .  ID - print id

   DESCRIPTION:
   This function print convergence rate.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

static INT PCR_IDAdmin, PCR_DispMode[32], PCR_nb[32], PCR_nComp[32], PCR_allComp[32], PCR_printed[32], PCR_nid[32];
static SHORT *PCR_ident[32];
static char PCR_compNames[32][MAX_VEC_COMP];
static const char *PCR_header[32];
static DOUBLE PCR_InitDefect[32][MAX_VEC_COMP], PCR_OldDefect[32][MAX_VEC_COMP];
static DOUBLE PCR_InitNorm[32],PCR_OldNorm[32];

static void PrintHeaderIff (INT i)
{
  INT j;

  if (PCR_header[i]==NULL)
    return;

  for (j=i+1; j<32; j++)
    if (PCR_printed[j])
      break;
  if (j<32)
  {
    UserWrite("\n");
    UserWrite(PCR_header[i]);
  }
}

INT PreparePCR (VECDATA_DESC *Vsym, INT DispMode, const char *text, INT *ID)
{
  INT i,j;

  /* get new ID */
  for (i=0; i<32; i++)
    if (((PCR_IDAdmin>>i)&1)==0)
    {
      PCR_IDAdmin |= (1<<i);
      *ID = i;
      break;
    }
  if (i==32)
  {
    PrintErrorMessage('E',"PreparePCR","no ID left");
    return (1);
  }

  /* init */
  PCR_nb[i]   = 0;
  PCR_DispMode[i] = DispMode;
  PCR_header[i] = text;
  for (j=i; j<32; j++) PCR_printed[j] = FALSE;

  /* print head line */
  if (text!=NULL && DispMode!=PCR_NO_DISPLAY)
  {
    UserWrite("\n");
    UserWrite(text);
  }

  /* store number and names of components */
  if (Vsym != NULL) {
    PCR_nComp[*ID] = VD_OFFSET(Vsym,NVECTYPES);
    if (PCR_nComp[*ID]>MAX_VEC_COMP) return (1);
    memcpy(PCR_compNames[*ID],VM_COMP_NAMEPTR(Vsym),MAX_VEC_COMP);
    PCR_nid[*ID] = VD_NID(Vsym);
    PCR_ident[*ID] = VD_IDENT_PTR(Vsym);
  }
  else if (*ID > 0) {
    PCR_nComp[*ID] = PCR_nComp[*ID-1];
    memcpy(PCR_compNames[*ID],PCR_compNames[*ID-1],MAX_VEC_COMP);
    PCR_nid[*ID] = PCR_nid[*ID-1];
    PCR_ident[*ID] = PCR_ident[*ID-1];
  }
  else {
    PCR_nComp[*ID] = MAX_VEC_COMP;
    memcpy(PCR_compNames[*ID],DEFAULT_NAMES,MAX_VEC_COMP);
    PCR_nid[*ID] = NO_IDENT;
  }
  PCR_allComp[*ID] = PCR_nComp[*ID];
  if (PCR_nid[*ID]!=NO_IDENT)
  {
    /* change compnames according to identification */
    for (j=0, i=0; i<PCR_nComp[*ID]; i++)
      if (PCR_ident[*ID][i]==i)
        PCR_compNames[*ID][j++] = PCR_compNames[*ID][i];

    PCR_nComp[*ID] = PCR_nid[*ID];
  }

  return (0);
}

/**************************************************************************/
/*
   PostPCR - stores the PCR values

   SYNOPSIS:
   INT PostPCR (INT ID, char *path);

   PARAMETERS:
   .  ID - print id
   .  path - envierement directory

   DESCRIPTION:
   This function stores the PCR values in the ug directory.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/*************************************************************************/

INT PostPCR (INT ID, char *path)
{
  INT i;
  char name[10];
  DOUBLE value,sum,defect;

  /* store values */
  if (path != NULL)
  {
    if (ChangeStructDir(path)==NULL) return (1);
    sum = 0.0;
    defect = 0.0;
    for (i=0; i<PCR_nComp[ID]; i++)
    {
      if (PCR_compNames[ID][i] == ' ') sprintf(name,"%c",(char)('a'+i));
      else sprintf(name,"%c",PCR_compNames[ID][i]);
      if (PCR_nb[ID]<=1) value = -1.0;
      else if (PCR_InitDefect[ID][i]==0.0) value = -2.0;
      else value = POW(PCR_OldDefect[ID][i]/PCR_InitDefect[ID][i],1.0/(PCR_nb[ID]-1));
      sum += value;
      defect += PCR_OldDefect[ID][i];
      if (SetStringValue(name,value)) return (1);
    }
    if (PCR_nComp[ID] > 0)
    {
      sum /= PCR_nComp[ID];
      if (SetStringValue("mean",sum)) return (1);
      defect /= PCR_nComp[ID];
      if (SetStringValue("defect",defect)) return (1);
    }
    if ((PCR_nComp[ID] > 1) && (PCR_InitNorm[ID] > 0) && (PCR_nb[ID] > 1))
      if (SetStringValue("norm",POW(PCR_OldNorm[ID]/PCR_InitNorm[ID],1.0/(PCR_nb[ID]-1)))) return (1);
    if (ChangeStructDir(":")==NULL) return (1);
  }

  if (ID>31 || ID<0 || (((PCR_IDAdmin>>ID)&1)==0)) return(1);
  PCR_IDAdmin &= ~(1<<(ID));
  return (0);
}

/*********************************************************************/
/*
   DoPCR - Print convergence rate routine

   SYNOPSIS:
   INT DoPCR (INT ID, VEC_SCALAR Defect, INT PrintMode);

   PARAMETERS:
   .  ID - print id
   .  Defect - array of DOUBLE
   .  PrintMode - mode

   DESCRIPTION:
   This function prints messages on the shell.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/**************************************************************************/

static INT NormIdentVS_of_VS (const VEC_SCALAR in, SHORT ncmp, SHORT nid, SHORT *ident, VEC_SCALAR out)
{
  INT i;

  if (nid!=NO_IDENT)
  {
    /* root of squared sum over identified components */
    DOUBLE sum;
    INT j,n=0;

    for (i=0; i<ncmp; i++)
      if (ident[i]==i)
      {
        sum = 0;
        for (j=0; j<ncmp; j++)
          if (ident[j]==i)
            sum += in[j]*in[j];
        out[n++] = sqrt(sum);
      }
    ASSERT(n==nid);
    return (0);
  }

  /* copy in --> out */
  for (i=0; i<ncmp; i++)
    out[i] = in[i];

  return (0);
}

INT DoPCR (INT ID, VEC_SCALAR InDefect, INT PrintMode)
{
  VEC_SCALAR Defect;
  DOUBLE d,s;
  INT i, j;

  /* check input */
  if (ID>31 || ID<0 || (((PCR_IDAdmin>>ID)&1)==0)) return(1);

  NormIdentVS_of_VS(InDefect,PCR_allComp[ID],PCR_nid[ID],PCR_ident[ID],Defect);

  /* calculate norm of defects */
  s = 0.0;
  for (j=0; j<PCR_nComp[ID]; j++)
  {
    d = Defect[j];
    s += d*d;
  }
  s = sqrt(s);

  switch (PrintMode)
  {
  case PCR_CRATE :
  case PCR_CRATE_SD :
    if (PCR_nb[ID]==0)
    {
      for (j=0; j<PCR_nComp[ID]; j++) PCR_InitDefect[ID][j] = Defect[j];
      PCR_InitNorm[ID] = s;
      if (PCR_DispMode[ID]==PCR_FULL_DISPLAY)
      {
        PCR_printed[ID] = TRUE;
        UserWriteF(" %-3d  %c: %-12.7e   %-12.7s\n",(int)PCR_nb[ID],PCR_compNames[ID][0],Defect[0],"---");
        for (i=1; i<PCR_nComp[ID]; i++)
          UserWriteF("      %c: %-12.7e   %-12.7s\n",PCR_compNames[ID][i],Defect[i],"---");
        if (PCR_nComp[ID]>1 && PrintMode==PCR_CRATE_SD)
          UserWriteF("   norm: %-12.7e   %-12.7s\n",s,"---");
        if (PCR_nComp[ID]>1) UserWrite("\n");
      }
    }
    else if (PCR_DispMode[ID]==PCR_FULL_DISPLAY)
    {
      PCR_printed[ID] = TRUE;
      PrintHeaderIff(ID);
      if (PCR_OldDefect[ID][0]!=0.0)
        UserWriteF(" %-3d  %c: %-12.7e   %-12.7e\n",(int)PCR_nb[ID],PCR_compNames[ID][0],Defect[0],Defect[0]/PCR_OldDefect[ID][0]);
      else
        UserWriteF(" %-3d  %c: %-12.7e   %-12.7s\n",(int)PCR_nb[ID],PCR_compNames[ID][0],Defect[0],"NaN");
      for (i=1; i<PCR_nComp[ID]; i++)
      {
        if (PCR_OldDefect[ID][i]!=0.0)
          UserWriteF("      %c: %-12.7e   %-12.7e\n",PCR_compNames[ID][i],Defect[i],Defect[i]/PCR_OldDefect[ID][i]);
        else
          UserWriteF("      %c: %-12.7e   %-12.7s\n",PCR_compNames[ID][i],Defect[i],"NaN");
      }
      if (PCR_nComp[ID]>1 && PrintMode==PCR_CRATE_SD)
        UserWriteF("   norm: %-12.7e   %-12.7e\n",s,s/PCR_OldNorm[ID]);
      if (PCR_nComp[ID]>1) UserWrite("\n");
    }
    for (j=0; j<PCR_nComp[ID]; j++) PCR_OldDefect[ID][j] = Defect[j];
    PCR_OldNorm[ID] = s;
    PCR_nb[ID]++;
    break;
  case PCR_AVERAGE :
  case PCR_AVERAGE_SD :
    if (PCR_nb[ID]<2) return (0);
    if (PCR_DispMode[ID]==PCR_NO_DISPLAY) break;
    PCR_printed[ID] = TRUE;
    PrintHeaderIff(ID);
    if (PCR_DispMode[ID]==PCR_FULL_DISPLAY) UserWrite("\n");
    if (PCR_InitDefect[ID][0]!=0.0)
      UserWriteF(" %-3d avg:  %c: %-12.7e   %-12.7e   %-12.7e\n",(int)(PCR_nb[ID]-1),PCR_compNames[ID][0],PCR_InitDefect[ID][0],Defect[0],POW(Defect[0]/PCR_InitDefect[ID][0],1.0/(PCR_nb[ID]-1)));
    else
      UserWriteF(" %-3d avg:  %c: %-12.7e   %-12.7e   %-12.7s\n",(int)(PCR_nb[ID]-1),PCR_compNames[ID][0],PCR_InitDefect[ID][0],Defect[0],"NaN");
    for (i=1; i<PCR_nComp[ID]; i++)
    {
      if (PCR_InitDefect[ID][i]!=0.0)
        UserWriteF("           %c: %-12.7e   %-12.7e   %-12.7e\n",PCR_compNames[ID][i],PCR_InitDefect[ID][i],Defect[i],POW(Defect[i]/PCR_InitDefect[ID][i],1.0/(PCR_nb[ID]-1)));
      else
        UserWriteF("           %c: %-12.7e   %-12.7e   %-12.7s\n",PCR_compNames[ID][i],PCR_InitDefect[ID][i],Defect[i],"NaN");
    }
    if (PCR_nComp[ID]>1 && PrintMode==PCR_AVERAGE_SD)
      UserWriteF("        norm: %-12.7e   %-12.7e   %-12.7e\n",PCR_InitNorm[ID],s,POW(s/PCR_InitNorm[ID],1.0/(PCR_nb[ID]-1)));
    UserWrite("\n");
    break;
  case PCR_INTERN :
  case PCR_INTERN_SD :
    PCR_nb[ID]++;
    s = 0.0;
    for (j=0; j<PCR_nComp[ID]; j++)
    {
      d = PCR_OldDefect[ID][j] = Defect[j];
      s += d*d;
    }
    s = sqrt(s);
    PCR_OldNorm[ID] = s;
    break;
  default :
    return (1);
  }

  return (0);
}
