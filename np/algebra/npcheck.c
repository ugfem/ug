// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      npcheck.c                                                     */
/*                                                                          */
/* Purpose:   check of numerical structures                                                     */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de			                                */
/*																			*/
/* History:   Juli 1 97 begin                                                                           */
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

#include <math.h>

#include "devices.h"
#include "compiler.h"
#include "gm.h"
#include "np.h"
#include "debug.h"
#include "general.h"

#include "npcheck.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


INT CheckSymmetryOfMatrix (GRID *theGrid, MATDATA_DESC *A)
{
  MATRIX *m,*mt;
  VECTOR *v,*w;
  DOUBLE *mptr,*mtptr;
  register SHORT i,j,rcomp,ccomp,*mcomp,*mtcomp,vtype,mtype;

  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    vtype = VTYPE(v);
    for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
      mt = MADJ(m);
      w = MDEST(m);
      mtype = MTP(vtype,VTYPE(w));
      rcomp = MD_ROWS_IN_MTYPE(A,mtype);
      if (rcomp == 0) continue;
      ccomp = MD_COLS_IN_MTYPE(A,mtype);
      if (ccomp == 0) continue;
      mcomp = MD_MCMPPTR_OF_MTYPE(A,mtype);
      mptr = MVALUEPTR(m,0);
      mtcomp = MD_MCMPPTR_OF_MTYPE(A,MTP(VTYPE(w),vtype));
      mtptr = MVALUEPTR(m,0);
      for (i=0; i<ccomp; i++)
        for (j=0; j<rcomp; j++)
          if (mptr[mcomp[i*rcomp+j]] != mtptr[mtcomp[j*ccomp+i]])
            return(1);
    }
  }

  return(0);
}

INT CheckVector(theGrid,theVector)
{
  INT nerr = 0;

  /* get format */

  /* check flags locally */

  return(nerr);
}

INT CheckVectors (GRID *theGrid)
{
  INT nerr = 0;
  VECTOR *theVector;

  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL; theVector=SUCCVC(theVector))
  {
    nerr += CheckVector(theGrid,theVector);
  }
  return(nerr);
}

INT CheckNP (MULTIGRID *theMG, INT argc, char **argv)
{
  MATDATA_DESC *A;
  INT level;
  char value[VALUELEN];

  if (ReadArgvChar("A",value,argc,argv) == 0) {
    A = GetMatDataDescByName(theMG,value);
    if (A != NULL)
      for (level=theMG->bottomLevel; level<=TOPLEVEL(theMG); level++)
        if (CheckSymmetryOfMatrix(GRID_ON_LEVEL(theMG,level),A))
          UserWriteF("matrix %s not symmetric on level %d\n",
                     ENVITEM_NAME(A),level);
  }

  for (level=theMG->bottomLevel; level<=TOPLEVEL(theMG); level++)
    if (CheckVectors(GRID_ON_LEVEL(theMG,level)))
      UserWriteF("ERROR: vector flags not correctly set");
  return(0);
}
