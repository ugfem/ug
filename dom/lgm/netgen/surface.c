// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      gg3d.c                                                        */
/*                                                                          */
/* Purpose:   interface for the 3d grid generator netgen                            */
/*                                                                          */
/* Author:    Christian Wieners                                             */
/*			  Institut fuer Computeranwendungen III                         */
/*			  Universitaet Stuttgart			                            */
/*			  Pfaffenwaldring 27				                            */
/*			  70569 Stuttgart, Germany			                            */
/*			  email: ug@ica3.uni-stuttgart.de		                        */
/*									                                        */
/* History:   18 March 96 begin, ug version 3.2                             */
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "compiler.h"
#include "devices.h"
#include "defaults.h"
#include "general.h"
#include "debug.h"
#include "lgm_domain.h"
#include "domain.h"


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

static INT nodeid;
static INT triangleid;
static INT left;
static INT right;
static double h_global;

static INT ntriangle;

static CoeffProcPtr LOCAL_H;

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/*****************************************************************************/


static INT AddPoint2Netgen (INT id, DOUBLE *global)
{

    #ifdef _NETGEN
  AddGeomPoint (nodeid,
                (double)global[0],(double)global[1],(double)global[2]);
    #endif

  return(0);
}


static LGM_SURFACE *theSurface;
int AddInnerNode2ug (double x, double y, double z)
{
  DOUBLE global[3],local[2];

  global[0] = x;
  global[1] = y;
  global[2] = z;

  /*printf("%s %f %f %f\n","outputpoint from netgen ",x,y,z);*/

  GetLocalKoord(theSurface,global,local);

  LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),0) = local[0];
  LGM_SURFDISC_LOCAL(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface)),1) = local[1];
  LGM_SURFDISC_NPOINT(LGM_SURFACE_DISC(theSurface))++;
  return(0);
}

int AddSurfaceTriangle2ug (int node0, int node1, int node2)
{
  INT Id[3];

  Id[0] = node0;
  Id[1] = node1;
  Id[2] = node2;

  /*printf("%s %d %d %d\n","outputtriangle from netgen ",Id[0],Id[1],Id[2]);*/

  LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)),0) = Id[0];
  LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)),1) = Id[1];
  LGM_SURFDISC_TRIANGLE(LGM_SURFACE_DISC(theSurface),LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface)),2) = Id[2];
  LGM_SURFDISC_NTRIANGLE(LGM_SURFACE_DISC(theSurface))++;

  return(0);
}


INT GenerateSurfaceGrid (LGM_SURFACE *aSurface, DOUBLE h, INT smooth,INT display)
{
  INT sid,i;
  char rulefilename[128];
  DOUBLE **x;

  theSurface = aSurface;
  ntriangle = 0;

  if (GetDefaultValue(DEFAULTSFILENAME,"netgentrianglerules",rulefilename))
    strcpy(rulefilename,rulefilename);

    #ifdef _NETGEN
  if (StartSurfaceNetgen(h,smooth,display)) return(1);
    #endif

  return(0);
}

int Get_Local_h(double *in, double *out)
{
  (*LOCAL_H)(in, out);
  return(0);
}

INT InitSurface(CoeffProcPtr Coeff)
{
  char rulefilename[128];
  if (GetDefaultValue(DEFAULTSFILENAME,"netgentrianglerules",rulefilename))
    strcpy(rulefilename,"triangle.rls");
  /*	LOCAL_H[0] = Coeff[coeff];*/
  LOCAL_H = Coeff;
    #ifdef _NETGEN
  InitSurfaceNetgen(rulefilename);
    #endif

}
