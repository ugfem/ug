// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      dummy-ppif.c                                                  */
/*                                                                          */
/* Purpose:   parallel processor interface                                  */
/*            Provides a portable interface to message passing MIMD         */
/*            architectures. PPIF is divided into three parts:              */
/*                                                                          */
/*            (1) Administration                                            */
/*            (2) Communication                                             */
/*            (3) Miscellaneous                                             */
/*                                                                          */
/*            The interface assumes that the parallel machine has           */
/*            the following properties:                                     */
/*                                                                          */
/*            (1) it is physically connected at least as a 2 or 3 dim. array*/
/*            (2) it has a fast virtual channel communication mechanism     */
/*            (3) it has an asynchronous communication mechanism            */
/*                                                                          */
/*            (Indigo) dummy module                                                     */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*            internet: bastian@iwr1.iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   17 Aug 1992, begin                                            */
/*            18 Feb 1993, Indigo version                                   */
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

/* standard C library */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "../ppif_general.h"

#include "compiler.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define MAXT        15      /* maximum number of downtree nodes max log2(P) */

/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
RCSID($Header$,PPIF_RCS_STRING)


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

typedef int *VChannelPtr;   /* dummy definition, any pointer type is ok     */

/* id's */
int me;                     /* my processor id                              */
int master;                 /* id of master processor                       */
int procs;                  /* number of processors in the network          */

/* 3D array structure */
int arrayid;                            /* compact format of position, 8 bits each      */
int MyX,MyY,MyZ;            /* 3D array coordinates                         */
int DimX,DimY,DimZ;         /* 3D array dimensions, may be 1 !              */
VChannelPtr nn[6];          /* nearest neighbors in 3D array                */

/* Tree structure */
int degree;                 /* degree of downtree nodes                     */
VChannelPtr uptree;         /* channel uptree                               */
VChannelPtr downtree[MAXT]; /* channels downtree (may be empty)             */
int slvcnt[MAXT];                       /* number of processors in subtree              */

/****************************************************************************/
/*                                                                          */
/* Function:  InitPPIF                                                      */
/*                                                                          */
/* Purpose:   initialize parallel processor interface                       */
/*            set exported variables, allocate tree communication structure */
/*                                                                          */
/* Input:     void                                                          */
/*                                                                          */
/* Output:    int 0:  ok                                                    */
/*            int!=0: error                                                 */
/*                                                                          */
/****************************************************************************/

int aid_to_pid (int x, int y, int z)
{
  if ((x<0)||(x>=DimX)) return(-1);
  if ((y<0)||(y>=DimY)) return(-1);
  if ((z<0)||(z>=DimZ)) return(-1);

  return((z*DimY+y)*DimX+x);
}

int pid_to_aid (int p)
{
  int x,y,z;

  if ((p<0)||(p>=procs)) return(-1);

  x = p%DimX;
  p = p/DimX;
  y = p%DimY;
  z = p/DimY;
  return((z<<16)|(y<<8)|x);
}


int InitPPIF (int *argcp, char ***argvp)
{
  int i;

  me = 0;
  master = 0;
  procs = 1;

  /* 3D array configuration */
  MyX = 0;
  MyY = 0;
  MyZ = 0;
  DimX = 1;
  DimY = 1;
  DimZ = 1;
  for (i=0; i<6; i++) nn[i] = NULL;

  /* tree configuration */
  degree = 0;
  uptree = NULL;
  for (i=0; i<MAXT; i++) downtree[i] = NULL;
  for (i=0; i<MAXT; i++) slvcnt[i] = 0;

  return(0);
}

void ExitPPIF (void)
{
  return;
}

/****************************************************************************/
/*                                                                          */
/* Tree oriented functions                                                  */
/*                                                                          */
/****************************************************************************/

int Broadcast (void *data, int size)
{
  return(0);
}

int Concentrate (void *data, int size)
{
  return(0);
}

int GetConcentrate (int slave, void *data, int size)
{
  return(0);
}

int Spread (int slave, void *data, int size)
{
  return(0);
}

int GetSpread (void *data, int size)
{
  return(0);
}

int Synchronize (void)
{
  return(0);
}


/****************************************************************************/
/*                                                                          */
/* Synchronous communication                                                */
/*                                                                          */
/****************************************************************************/

VChannelPtr ConnSync (int p, int id)
{
  return(0);
}

int DiscSync (VChannelPtr vc)
{
  return(0);
}

int SendSync (VChannelPtr vc, void *data, int size)
{
  return(size);
}

int RecvSync (VChannelPtr vc, void *data, int size)
{
  return(size);
}


/****************************************************************************/
/*                                                                          */
/* Asynchronous communication                                               */
/*                                                                          */
/****************************************************************************/

int SendASync (VChannelPtr vc, void *data, int size)
{
  return(0);
}

int RecvASync (VChannelPtr vc, void *data, int size)
{
  return(0);
}

int InfoASend (VChannelPtr vc)
{
  return(0);
}

int InfoARecv (VChannelPtr vc)
{
  return(0);
}

VChannelPtr ConnASync (int p, int id)
{
  return((VChannelPtr) 4);
}

int InfoAConn (VChannelPtr vc)
{
  return(0);
}


int DiscASync (VChannelPtr vc)
{
  return(0);
}

int InfoADisc (VChannelPtr vc)
{
  return(0);
}


/****************************************************************************/
/*                                                                          */
/* Random communication                                                         */
/*                                                                          */
/****************************************************************************/

int SendMail (int destId, int reqId, void *data, int size)
{
  return(0);
}

int GetMail (int *sourceId, int *reqId, void *data, int *size)
{
  return(0);
}


/****************************************************************************/
/*                                                                          */
/* Miscellaneous                                                                        */
/*                                                                          */
/****************************************************************************/

int UsedSpace (void)
{
  return(0);
}

void PrintHostMessage (char *s)
{
  printf("%s",s);
}

double CurrentTime (void)
{
  return(((double)(clock())/((double)CLOCKS_PER_SEC)));
}

int Distance (int p, int q)
{
  return(0);
}
