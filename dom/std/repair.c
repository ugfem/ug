// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  repair.c	                                                                                                */
/*																			*/
/* Purpose:   repair tool                                                                               */
/*																			*/
/* Author:	  Christian Wieners                                                                     */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   Jun 28 for Marc                                                                   */
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

/* standard C library */
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

/* low modules */
#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"
#include "bio.h"
#include "misc.h"
#include "fileopen.h"
#include "defaults.h"
#include "general.h"
#include "debug.h"
#include "evm.h"

/* dev modules */
#include "ugdevices.h"

/* domain module */
#include "std_domain.h"
#include "domain.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define SMALL_DIFF   SMALL_C*100
#define RESOLUTION   100

#define DEFAULTDOMMEMORY 50000

#define OPTIONLEN 32

#define V2_LINCOMB(a,A,b,B,C)              {(C)[0] = (a)*(A)[0] + (b)*(B)[0];\
                                            (C)[1] = (a)*(A)[1] + (b)*(B)[1];}

#define V2_EUKLIDNORM_OF_DIFF(A,B,b)    (b) = sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])+((A)[1]-(B)[1])*((A)[1]-(B)[1])));

#define V3_EUKLIDNORM_OF_DIFF(A,B,b)    (b) = (sqrt((double)(((A)[0]-(B)[0])*((A)[0]-(B)[0])+((A)[1]-(B)[1])*((A)[1]-(B)[1])+((A)[2]-(B)[2])*((A)[2]-(B)[2]))));

#ifdef __TWODIM__
#define V_DIM_EUKLIDNORM_OF_DIFF(A,B,b) V2_EUKLIDNORM_OF_DIFF(A,B,b)
#else
#define V_DIM_EUKLIDNORM_OF_DIFF(A,B,b) V3_EUKLIDNORM_OF_DIFF(A,B,b)
#endif

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

static INT theProblemDirID;             /* env type for Problem dir                     */
static INT theBdryCondVarID;            /* env type for Problem vars			*/

static INT theDomainDirID;                      /* env type for Domain dir				*/
static INT theBdrySegVarID;             /* env type for bdry segment vars		*/
static INT theLinSegVarID;                  /* env type for linear segment vars		*/

static INT theBVPDirID;                         /* env type for BVP dir					*/

static STD_BVP *currBVP;

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

#ifdef __THREEDIM__

typedef struct {
  INT corners_of_elem;                                          /* number of corners                    */
  INT sides_of_elem;                                                    /* number of sides                      */
  INT corners_of_side[MAX_SIDES_OF_ELEM];       /* nb. of corners for each side */
  INT corner_of_side[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE];

} REFERENCE_ELEMENT;

static REFERENCE_ELEMENT Tetrahedron = {
  4,                                                                                    /* tag							*/
  4,                                                                                    /* number of sides				*/
  {3,3,3,3,-1,-1},                                                      /* corners for each side		*/
  {{0,2,1,-1},{1,2,3,-1},{0,3,2,-1},{0,1,3,-1}}
};

static REFERENCE_ELEMENT Pyramid = {
  5,                                                                                    /* tag							*/
  5,                                                                                    /* number of sides				*/
  {4,3,3,3,3,-1},                                                       /* corners for each side		*/
  {{0,3,2,1},{0,1,4,-1},{1,2,4,-1},                     /* number of corner j of side i */
   {2,3,4,-1},{3,0,4,-1}}
};

static REFERENCE_ELEMENT Prism = {
  6,                                                                                    /* tag							*/
  5,                                                                                    /* number of sides				*/
  {3,4,4,4,3,-1},                                                       /* corners for each side		*/
  {{0,2,1,-1},{0,1,4,3},{1,2,5,4},                      /* number of corner j of side i */
   {2,0,3,5},{3,4,5,-1}}
};

static REFERENCE_ELEMENT Hexahedron = {
  8,                                                                                    /* tag							*/
  6,                                                                                    /* number of sides				*/
  {4,4,4,4,4,4},                                                        /* corners for each side		*/
  {{0,3,2,1},{0,1,5,4},{1,2,6,5},                       /* number of corner j of side i */
   {2,3,7,6},{3,0,4,7},{4,5,6,7}}
};


static INT CheckOnSide(INT *corner, INT n, INT **ids, INT *flag)
{
  INT rv = 0;
  INT i,j,k,m;

  for (j=0; j<6; j++)
  {
    INT side[4];

    for (k=0; k<4; k++)
      side[k] = corner[Hexahedron.corner_of_side[j][k]];

    flag[j] = 0;

    for (i=0; i<n; i++)
    {
      m = 0;
      for (k=0; k<4; k++) {
        if (side[k] == ids[i][Prism.corner_of_side[0][0]]) m++;
        else if (side[k] == ids[i][Prism.corner_of_side[0][1]]) m++;
        else if (side[k] == ids[i][Prism.corner_of_side[0][2]]) m++;
      }
      if (m == 3) {
        rv = 1;
        for (k=0; k<4; k++) {
          if (side[k] == ids[i][Prism.corner_of_side[0][0]])
            continue;
          else if (side[k] == ids[i][Prism.corner_of_side[0][1]])
            continue;
          else if (side[k] == ids[i][Prism.corner_of_side[0][2]])
            continue;
          break;
        }
        if (k == 0) flag[j] = 2;
        else if (k == 1) flag[j] = 1;
        else if (k == 2) flag[j] = 2;
        else if (k == 3) flag[j] = 1;


        if (0) {
          UserWriteF("\nside[%d]  ",j);
          for(k=0; k<4; k++)
            UserWriteF("_%05d",side[k]);
          UserWriteF("\n");
        }

      }
      m = 0;
      for (k=0; k<4; k++) {
        if (side[k] == ids[i][Prism.corner_of_side[4][0]]) m++;
        else if (side[k] == ids[i][Prism.corner_of_side[4][1]]) m++;
        else if (side[k] == ids[i][Prism.corner_of_side[4][2]]) m++;
      }
      if (m == 3) {
        rv = 1;
        for (k=0; k<4; k++) {
          if (side[k] == ids[i][Prism.corner_of_side[4][0]])
            continue;
          else if (side[k] == ids[i][Prism.corner_of_side[4][1]])
            continue;
          else if (side[k] == ids[i][Prism.corner_of_side[4][2]])
            continue;
          break;
        }
        if (k == 0) assert(flag[j] == 2);
        else if (k == 1) assert(flag[j] == 1);
        else if (k == 2) assert(flag[j] == 2);
        else if (k == 3) assert(flag[j] == 1);
      }
    }
  }

  return(rv);
}

static INT CheckOrientation (INT i, INT n,
                             DOUBLE *x, DOUBLE *y, DOUBLE *z, DOUBLE *w)
{
  DOUBLE_VECTOR diff[3],rot;
  DOUBLE det;
  INT j;

  V3_SUBTRACT(y,x,diff[0]);
  V3_SUBTRACT(z,x,diff[1]);
  V3_SUBTRACT(w,x,diff[2]);
  V3_VECTOR_PRODUCT(diff[0],diff[1],rot);
  V3_SCALAR_PRODUCT(rot,diff[2],det);

  if (det < 0.0) {
    UserWriteF(" ID(Elem)=%d n=%d  det %f\n",i,n,det);
    for (j=0; j<3; j++)
      UserWriteF("        diff[%d]=%5.2f %5.2f %5.2f\n",
                 j,diff[j][0],diff[j][1],diff[j][2]);
    UserWriteF("\n");

    return(1);
  }
  return(0);
}
#endif

INT RepairMesh (HEAP *Heap, INT MarkKey, MESH *Mesh)
{
  INT i,sd;
  DOUBLE **pos = Mesh->Position;
  DOUBLE **p;
  INT c[8];

#ifdef __THREEDIM__

  Mesh->Position = (DOUBLE **)
                   GetTmpMem(Heap,Mesh->nInnP*2*sizeof(DOUBLE *),MarkKey);
  for (i=0; i<Mesh->nInnP; i++)
    Mesh->Position[i] = pos[i];

  p = (DOUBLE **)
      GetTmpMem(Heap,(Mesh->nInnP+Mesh->nBndP)*sizeof(DOUBLE *),MarkKey);

  for (i=0; i<Mesh->nBndP; i++) {
    p[i] = (DOUBLE *)
           GetTmpMem(Heap,3*sizeof(DOUBLE),MarkKey);
    BNDP_Global(Mesh->theBndPs[i],p[i]);
  }
  for (i=Mesh->nBndP; i<Mesh->nBndP+Mesh->nInnP; i++)
    p[i] = pos[i-Mesh->nBndP];

  for (sd=1; sd<=Mesh->nSubDomains; sd++)
  {
    INT nElem = Mesh->nElements[sd];
    INT *corners = Mesh->Element_corners[sd];
    INT **corner = Mesh->Element_corner_ids[sd];
    INT j,k,nPrism;


    if (0)
      for (i=0; i<nElem; i++) {
        if (corners[i] == 6) {
          CheckOrientation(i,6,
                           p[corner[i][0]],p[corner[i][1]],
                           p[corner[i][2]],p[corner[i][3]]);

          /*
             UserWriteF(" ID(Elem)=%d\n",i);
             for (j=0; j<6; j++)
             UserWriteF("        pos[%d]=%5.2f %5.2f %5.2f\n",
                                   j,p[corner[i][j]][0],
                                   p[corner[i][j]][1],
                                   p[corner[i][j]][2]);
           */

        }
        if (corners[i] == 8) {
          if (CheckOrientation(i,8,
                               p[corner[i][0]],p[corner[i][1]],
                               p[corner[i][2]],p[corner[i][4]])) {

            for (j=0; j<8; j++)
              c[j] = corner[i][j];
            for (j=0; j<4; j++) {
              corner[i][j] = c[j+4];
              corner[i][j+4] = c[j];
            }
          }
          CheckOrientation(i,8,
                           p[corner[i][0]],p[corner[i][1]],
                           p[corner[i][2]],p[corner[i][4]]);


          /*
             UserWriteF(" ID(Elem)=%d\n",i);
             for (j=0; j<8; j++)
             UserWriteF("        pos[%d]=%5.2f %5.2f %5.2f\n",
                                   j,p[corner[i][j]][0],
                                   p[corner[i][j]][1],
                                   p[corner[i][j]][2]);
           */

        }
      }

    nPrism = 0;
    for (i=0; i<nElem; i++)
      if (corners[i] == 6) nPrism++;

    Mesh->Element_corners[sd] = (INT *)
                                GetTmpMem(Heap,(nElem+11*(nElem-nPrism))*sizeof(INT),MarkKey);
    Mesh->Element_corner_ids[sd] = (INT **)
                                   GetTmpMem(Heap,(nElem+11*(nElem-nPrism))*sizeof(INT*),MarkKey);

    nPrism = 0;
    for (i=0; i<nElem; i++)
      if (corners[i] == 6) {
        Mesh->Element_corners[sd][nPrism] = 6;
        Mesh->Element_corner_ids[sd][nPrism] = corner[i];
        /*
           Mesh->Element_corner_ids[sd][nPrism] = (INT *)
            GetTmpMem(Heap,6*sizeof(INT),MarkKey);
           for (j=0; j<6; j++)
            Mesh->Element_corner_ids[sd][nPrism][j] = corner[i][j];
         */
        nPrism++;
      }
    Mesh->nElements[sd] = nPrism;

    for (i=0; i<nElem; i++)
    {
      INT flag[6];

      if (corners[i] == 6) continue;

      if (0) {

        UserWriteF("repair element[%d]\n",i);
        for (j=0; j<6; j++)
        {
          UserWriteF("corner[%d]  ",j);
          for(k=0; k<4; k++)
          {
            UserWriteF("_%05d",corner[i][Hexahedron.corner_of_side[j][k]]);
          }
          UserWriteF("\n");
        }
        UserWriteF("\n");
      }


      if (CheckOnSide(corner[i],nPrism,
                      Mesh->Element_corner_ids[sd],flag))
      {
        DOUBLE_VECTOR global;
        DOUBLE scale = 0.125;

        V_DIM_CLEAR(global);

        for (j=0; j<8; j++) {
          V_DIM_LINCOMB(1.0,global,scale,pos[corner[i][j]],global);
        }
        Mesh->Position[Mesh->nInnP] = (DOUBLE *)
                                      GetTmpMem(Heap,3*sizeof(DOUBLE),MarkKey);
        V_DIM_COPY(global,Mesh->Position[Mesh->nInnP]);


        if (0) {
          UserWriteF(" el[%i] midPoint %8.4f %8.4f %8.4f\n",
                     i,global[0],global[1],global[2]);
          for (j=0; j<8; j++)
            UserWriteF("   corner[%d] %8.4f %8.4f %8.4f\n",
                       j,pos[corner[i][j]][0],
                       pos[corner[i][j]][1],
                       pos[corner[i][j]][2]);
        }


        for (j=0; j<6; j++) {
          if (flag[j] == 0) {
            Mesh->Element_corners[sd][Mesh->nElements[sd]] = 5;
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]] =
              (INT *) GetTmpMem(Heap,5*sizeof(INT),MarkKey);
            for (k=0; k<4; k++)
              Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][Pyramid.corner_of_side[0][k]]=
                corner[i][Hexahedron.corner_of_side[j][k]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][4] =
              Mesh->nInnP;

            if(0)
            {
              INT j1,k1;


              UserWriteF("repaired element[%d]\n",Mesh->nElements[sd]);
              for (j1=0; j1<5; j1++)
              {
                UserWriteF("corner[%d]  ",j1);
                for(k1=0; k1<Pyramid.corners_of_side[j1]; k1++)
                {
                  UserWriteF("_%05d",
                             Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][Pyramid.corner_of_side[j1][k1]]);
                }
                UserWriteF("\n");
              }
              UserWriteF("\n");

            }


            Mesh->nElements[sd]++;
          }
          else if (flag[j] == 1) {
            Mesh->Element_corners[sd][Mesh->nElements[sd]] = 4;
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]] =
              (INT *) GetTmpMem(Heap,4*sizeof(INT),MarkKey);
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][0]=
              corner[i][Hexahedron.corner_of_side[j][0]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][1]=
              corner[i][Hexahedron.corner_of_side[j][2]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][2]=
              corner[i][Hexahedron.corner_of_side[j][1]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][3] =
              Mesh->nInnP;



            if(0)
            {
              INT j1,k1;


              UserWriteF("A repaired element[%d] flag %d\n",Mesh->nElements[sd],flag[j]);
              for (j1=0; j1<4; j1++)
              {
                UserWriteF("corner[%d]  ",j1);
                for(k1=0; k1<Tetrahedron.corners_of_side[j1]; k1++)
                {
                  UserWriteF("_%05d",
                             Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][Tetrahedron.corner_of_side[j1][k1]]);
                }
                UserWriteF("\n");
              }
              UserWriteF("\n");

            }









            Mesh->nElements[sd]++;
            Mesh->Element_corners[sd][Mesh->nElements[sd]] = 4;
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]] =
              (INT *) GetTmpMem(Heap,4*sizeof(INT),MarkKey);
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][0]=
              corner[i][Hexahedron.corner_of_side[j][0]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][2]=
              corner[i][Hexahedron.corner_of_side[j][2]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][1]=
              corner[i][Hexahedron.corner_of_side[j][3]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][3] =
              Mesh->nInnP;

            if(0)
            {
              INT j1,k1;


              UserWriteF("B repaired element[%d] flag %d\n",Mesh->nElements[sd],flag[j]);
              for (j1=0; j1<4; j1++)
              {
                UserWriteF("corner[%d]  ",j1);
                for(k1=0; k1<Tetrahedron.corners_of_side[j1]; k1++)
                {
                  UserWriteF("_%05d",
                             Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][Tetrahedron.corner_of_side[j1][k1]]);
                }
                UserWriteF("\n");
              }
              UserWriteF("\n");

            }



            Mesh->nElements[sd]++;
          }
          else if (flag[j] == 2) {
            Mesh->Element_corners[sd][Mesh->nElements[sd]] = 4;
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]] =
              (INT *) GetTmpMem(Heap,4*sizeof(INT),MarkKey);
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][0]=
              corner[i][Hexahedron.corner_of_side[j][0]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][2]=
              corner[i][Hexahedron.corner_of_side[j][1]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][1]=
              corner[i][Hexahedron.corner_of_side[j][3]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][3] =
              Mesh->nInnP;

            if(0)
            {
              INT j1,k1;


              UserWriteF("C repaired element[%d] flag %d\n",Mesh->nElements[sd],flag[j]);
              for (j1=0; j1<4; j1++)
              {
                UserWriteF("corner[%d]  ",j1);
                for(k1=0; k1<Tetrahedron.corners_of_side[j1]; k1++)
                {
                  UserWriteF("_%05d",
                             Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][Tetrahedron.corner_of_side[j1][k1]]);
                }
                UserWriteF("\n");
              }
              UserWriteF("\n");

            }



            Mesh->nElements[sd]++;
            Mesh->Element_corners[sd][Mesh->nElements[sd]] = 4;
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]] =
              (INT *) GetTmpMem(Heap,4*sizeof(INT),MarkKey);
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][0]=
              corner[i][Hexahedron.corner_of_side[j][1]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][2]=
              corner[i][Hexahedron.corner_of_side[j][2]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][1]=
              corner[i][Hexahedron.corner_of_side[j][3]];
            Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][3] =
              Mesh->nInnP;

            if(0)
            {
              INT j1,k1;


              UserWriteF("D repaired element[%d] flag %d\n",Mesh->nElements[sd],flag[j]);
              for (j1=0; j1<4; j1++)
              {
                UserWriteF("corner[%d]  ",j1);
                for(k1=0; k1<Tetrahedron.corners_of_side[j1]; k1++)
                {
                  UserWriteF("_%05d",
                             Mesh->Element_corner_ids[sd][Mesh->nElements[sd]][Tetrahedron.corner_of_side[j1][k1]]);
                }
                UserWriteF("\n");
              }
              UserWriteF("\n");

            }



            Mesh->nElements[sd]++;
          }
        }
        Mesh->nInnP++;
      }
      else {
        Mesh->Element_corners[sd][Mesh->nElements[sd]] = 8;
        Mesh->Element_corner_ids[sd][Mesh->nElements[sd]] = corner[i];
        Mesh->nElements[sd]++;
      }
    }
  }
#endif

  return(0);
}
