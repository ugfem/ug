// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      mpi-ppif.c                                                    */
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
/*            MPI Implementation                                            */
/*                                                                          */
/* Author:    Jens Boenisch                                                 */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 3a                                                */
/*            7000 Stuttgart 80                                             */
/*            internet: boenisch@rus.uni-stuttgart.de                       */
/*                                                                          */
/* History:   17 Aug 1992, begin                                            */
/*            18 Feb 1993, Indigo version                                   */
/*            02 Jun 1993, Paragon OSF/1 version                            */
/*            14 Sep 1995, MPI version
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
#include <sys/types.h>
#include <curses.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>



/* MPI library */
/*
   #include <mpp/mpi.h>
 */
#include <mpi.h>

#include "../ppif_general.h"

/*#include "compiler.h"*/

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define RAND_MSG_SIZE   128     /* max size of random messages              */
#define MAXT            15      /* maximum number of downtree nodes max     */
                                /* log2(P)                                  */
#define MAXVCHAN        256     /* maximum number of channels               */
#define PTYPE_ANY       -1L     /* process type: any process                */

#define ID_ARRAY        100L    /* channel id: array                        */
#define ID_TREE         101L    /* channel id: tree                         */
#define ID_GLOBAL       102L    /* channel id: global                       */
#define ID_MAIL         103L    /* channel id: mail                         */

#define ABS(i)          (((i)<0) ? (-(i)) : (i))

#define XPOS(aid)       (aid&0xFF)              /* xpos from compact form   */
#define YPOS(aid)       ((aid&0xFF00)>>8)       /* ypos from compact form   */
#define ZPOS(aid)       ((aid&0xFF0000)>>16)    /* zpos from compact form   */

#define PPIF_SUCCESS    0L      /* Return value for success                 */
#define PPIF_FAILURE    1L      /* Return value for failure            `    */

/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

enum directions {north, east, south, west, up, down};

typedef struct {
  int used;
  long p;
  long chanid;
} VChannel;

typedef VChannel *VChannelPtr;
typedef int      *msgid;


/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
/*RCSID("$Header$",PPIF_RCS_STRING)*/


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/* id's */
int me;                         /* my processor id                          */
int master;                     /* id of master processor                   */
int procs;                      /* number of processors in the network      */

/* 3D array structure */
int arrayid;                    /* compact format of position, 8 bits each  */
int MyX, MyY, MyZ;              /* 3D array coordinates                     */
int DimX, DimY, DimZ;           /* 3D array dimensions, may be 1 !          */
VChannelPtr nn[6];              /* nearest neighbors in 3D array            */

/* Tree structure */
int degree;                     /* degree of downtree nodes                 */
VChannelPtr uptree;             /* channel uptree                           */
VChannelPtr downtree[MAXT];     /* channels downtree (may be empty)         */
int slvcnt[MAXT];               /* number of processors in subtree          */

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static VChannel *VChan;
static int vc_count=0;          /* number of used VChan                     */
static int vc_free=0;           /* rotating pointer to find free VChan      */

static char buf[100];           /* string buffer for debugging              */


/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

int SendSync (VChannelPtr vc, void *data, int size);
int RecvSync (VChannelPtr vc, void *data, int size);

/****************************************************************************/
/*                                                                          */
/* routines for handling virtual channels                                                           */
/*                                                                          */
/****************************************************************************/

int InitVChan (void)
{
  int i;

  if (VChan = (VChannel *) malloc(sizeof(VChannel)*MAXVCHAN) )
  {
    for (i=0; i<MAXVCHAN; i++) VChan[i].used = FALSE;

    vc_count=0;
    vc_free=0;
  }

  return ((int) VChan);
}


VChannelPtr NewVChan (int p, int id)

{
  int i;
  VChannelPtr myChan;

  if (vc_count < MAXVCHAN)
  {
    while (VChan[vc_free].used)
    {
      if (++vc_free == MAXVCHAN) vc_free=0;
    }

    myChan = &VChan[vc_free];

    myChan->used   = TRUE;
    myChan->p      = p;
    myChan->chanid = id;

    vc_count += 1;
    vc_free = (vc_free + 1) % MAXVCHAN;

    return (myChan);
  }
  else
  {
    printf ("%d: no more VChannels in NewVChan(), dest=%d, id=%d\n",me,p,id);

    return (NULL);
  }
}


void DeleteVChan (VChannelPtr myChan)

{
  myChan->used = FALSE;

  vc_count -= 1;
}


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
  if ((x<0)||(x>=DimX)) return (-1);
  if ((y<0)||(y>=DimY)) return (-1);
  if ((z<0)||(z>=DimZ)) return (-1);

  return ( (z*DimY+y)*DimX+x);
}

int pid_to_aid (int p)

{
  int x, y, z;

  if ((p<0)||(p>=procs)) return (-1);

  x = p%DimX;
  p = p/DimX;
  y = p%DimY;
  z = p/DimY;

  return ((z<<16)|(y<<8)|x);
}

/*
   Factor N into two integers that are as close together as possible
 */

void Factor (int N, int *pn, int *pm)

{
  int n = ceil (sqrt ((double) N));
  int m = floor (sqrt ((double) N));

  while (n*m != N)
  {
    if (n*m < N) n++;
    else m--;
  }

  *pn = n; *pm = m;
}


int InitPPIF (int *argcp, char ***argvp)
{
  int i, succ, sonr, sonl, aid;
  MPI_Status status;

  MPI_Init (argcp, argvp);

  MPI_Comm_rank (MPI_COMM_WORLD, &me);
  MPI_Comm_size (MPI_COMM_WORLD, &procs);


  master = 0;

  if (!InitVChan()) printf(" %4d: Couldn't get VChannel memory!\n", me);

  DimZ = 1;
  Factor(procs, &DimX, &DimY);

  if (me==master) printf("DimX=%d, DimY=%d, DimZ=%d\n", DimX, DimY, DimZ);

  aid = pid_to_aid(me);
  MyX = XPOS(aid);
  MyY = YPOS(aid);
  MyZ = ZPOS(aid);
  for (i=0; i<6; i++) nn[i] = NULL;

  /* alloc east-west */
  if (DimX>1)
  {
    if (MyX>0)
    {
      nn[west] = NewVChan(me-1,ID_ARRAY);
      if (nn[west]==NULL) return(1);
    }
    if (MyX<DimX-1)
    {
      nn[east] =NewVChan(me+1,ID_ARRAY);
      if (nn[east]==NULL) return(1);
    }
  }

  /* alloc north-south */
  if (DimY>1)
  {
    if (MyY>0)
    {
      nn[south] = NewVChan(me-DimX,ID_ARRAY);
      if (nn[south]==NULL) return(1);
    }
    if (MyY<DimY-1)
    {
      nn[north] =NewVChan(me+DimX,ID_ARRAY);
      if (nn[north]==NULL) return(1);
    }
  }

  /* alloc up-down */
  if (DimZ>1)
  {
    if (MyZ>0)
    {
      nn[down] = NewVChan(me-DimX*DimY,ID_ARRAY);
      if (nn[down]==NULL) return(1);
    }
    if (MyZ<DimZ-1)
    {
      nn[up] = NewVChan(me+DimX*DimY,ID_ARRAY);
      if (nn[up]==NULL) return(1);
    }
  }

  /* tree configuration */

  degree = 0;
  sonl = 2*me + 1;
  sonr = 2*me + 2;

  if (sonl<procs)
  {
    degree++;
    downtree[0] = NewVChan(sonl,ID_TREE);
  }
  else
  {
    downtree[0] = NULL;
  }

  if (sonr<procs)
  {
    degree++;
    downtree[1] = NewVChan(sonr,ID_TREE);
  }
  else
  {
    downtree[1] = NULL;
  }

  if (me>0)
  {
    uptree = NewVChan((me-1)/2,ID_TREE);
  }
  else
  {
    uptree = NULL;
  }

  succ=1;
  for(i=0; i<degree; i++)
  {
    MPI_Recv ((void *) &(slvcnt[i]), (long) sizeof(int), MPI_BYTE, downtree[i]->p, ID_TREE, MPI_COMM_WORLD, &status);
    succ += slvcnt[i];
  }
  if (me>0)
  {
    MPI_Send ((void *) &succ, (long) sizeof(int), MPI_BYTE, (long)(me-1)/2, ID_TREE, MPI_COMM_WORLD);
  }

  return (PPIF_SUCCESS);
}

void ExitPPIF (void)

{
  MPI_Finalize ();
}

/****************************************************************************/
/*                                                                          */
/* Tree oriented functions                                                  */
/*                                                                          */
/****************************************************************************/

int Broadcast (void *data, int size)

{
  if (MPI_SUCCESS != MPI_Bcast (data, size, MPI_BYTE, master, MPI_COMM_WORLD) )
    return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int Concentrate (void *data, int size)

{
  if (me != master)
    if (SendSync (uptree, data, size) < 0) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int GetConcentrate (int slave, void *data, int size)

{
  VChannel *chan;

  if (slave < degree)
    if (RecvSync (downtree[slave], data, size) < 0) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int Spread (int slave, void *data, int size)

{
  if (slave < degree)
    if (SendSync (downtree[slave], data, size) < 0) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int GetSpread (void *data, int size)

{
  if (me!=master)
    if (RecvSync (uptree, data, size) < 0) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int Synchronize (void)

{
  if (MPI_SUCCESS != MPI_Barrier (MPI_COMM_WORLD) ) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}


/****************************************************************************/
/*                                                                          */
/* Synchronous communication                                                */
/*                                                                          */
/****************************************************************************/

VChannelPtr ConnSync (int p, int id)

{
  return (NewVChan (p, id) );
}

int DiscSync (VChannelPtr vc)

{
  DeleteVChan (vc);

  return (0);
}

int SendSync (VChannelPtr vc, void *data, int size)

{
  if (MPI_SUCCESS == MPI_Ssend (data, size, MPI_BYTE, vc->p, vc->chanid, MPI_COMM_WORLD) )
    return (size);
  else
    return (-1);
}

int RecvSync (VChannelPtr vc, void *data, int size)

{
  int count = -1;
  MPI_Status status;

  if (MPI_SUCCESS == MPI_Recv (data, size, MPI_BYTE, vc->p, vc->chanid, MPI_COMM_WORLD, &status) )
    MPI_Get_count (&status, MPI_BYTE, &count);

  return (count);
}


/****************************************************************************/
/*                                                                          */
/* Asynchronous communication                                               */
/*                                                                          */
/****************************************************************************/

VChannelPtr ConnASync (int p, int id)

{
  return (NewVChan (p,id) );
}

int InfoAConn (VChannelPtr vc)

{
  return (vc ? 1 : -1);
}


int DiscASync (VChannelPtr vc)

{
  DeleteVChan (vc);
  return (PPIF_SUCCESS);
}

int InfoADisc (VChannelPtr vc)

{
  return (TRUE);
}

#define REQUEST_HEAP

msgid SendASync (VChannelPtr vc, void *data, int size, int *error)

{
#  ifdef REQUEST_HEAP
  MPI_Request *req;

  if (req = malloc (sizeof (MPI_Request) ) )
  {
    if (MPI_SUCCESS == MPI_Isend (data, size, MPI_BYTE, vc->p, vc->chanid, MPI_COMM_WORLD, req) )
    {
      *error = FALSE;
      return ((msgid) req);
    }
  }

#  else
  MPI_Request Req;
  if (MPI_SUCCESS == MPI_Isend (data, size, MPI_BYTE, vc->p, vc->chanid, MPI_COMM_WORLD, &Req) )
  {
    *error = FALSE;
    return ((msgid) Req);
  }
#  endif

  *error = TRUE;
  return NULL;
}


msgid RecvASync (VChannelPtr vc, void *data, int size, int *error)

{
#  ifdef REQUEST_HEAP
  MPI_Request *req;

  if (req = malloc (sizeof (MPI_Request) ) )
  {
    if (MPI_SUCCESS == MPI_Irecv (data, size, MPI_BYTE, vc->p, vc->chanid, MPI_COMM_WORLD, req) )
    {
      *error = FALSE;
      return ((msgid) req);
    }
  }

#  else
  MPI_Request Req;

  if (MPI_SUCCESS == MPI_Irecv (data, size, MPI_BYTE, vc->p, vc->chanid, MPI_COMM_WORLD, &Req) )
  {
    *error = FALSE;
    return ((msgid) Req);
  }

#  endif

  *error = TRUE;
  return (NULL);
}


int InfoASend (VChannelPtr vc, msgid m)

{
  MPI_Status status;
  int complete;

#  ifdef REQUEST_HEAP
  if (m)
  {
    if (MPI_SUCCESS == MPI_Test ((MPI_Request *) m, &complete, &status) )
    {
      if (complete) free (m);

      return (complete);        /* complete is TRUE for completed send, FALSE otherwise */
    }
  }

#  else
  MPI_Request Req = (MPI_Request) m;

  if (MPI_SUCCESS == MPI_Test (&Req, &complete, &status) )
  {
    return (complete);          /* complete is TRUE for completed send, FALSE otherwise */
  }

#  endif

  return (-1);          /* return -1 for FAILURE */
}


int InfoARecv (VChannelPtr vc, msgid m)
{
  MPI_Status status;
  int complete;

#  ifdef REQUEST_HEAP
  if (m)
  {
    if (MPI_SUCCESS == MPI_Test ((MPI_Request *) m, &complete, &status) )
    {
      if (complete) free (m);

      return (complete);        /* complete is TRUE for completed receive, FALSE otherwise */
    }
  }

#  else
  MPI_Request Req = (MPI_Request) m;

  if (MPI_SUCCESS == MPI_Test (&Req, &complete, &status) )
  {
    return (complete);          /* complete is TRUE for completed receive, FALSE otherwise */
  }

#  endif

  return (-1);          /* return -1 for FAILURE */
}


/****************************************************************************/
/*                                                                          */
/* Random communication                                                         */
/*                                                                          */
/****************************************************************************/

int SendMail (int destId, int reqId, void *data, int size)

{
  if (MPI_SUCCESS == MPI_Ssend (data, size, MPI_BYTE, destId, ID_MAIL, MPI_COMM_WORLD) )
    return (PPIF_SUCCESS);

  return (PPIF_FAILURE);
}

int GetMail (int *sourceId, int *reqId, void *data, int *size)

{
  MPI_Status status;
  int flag;

  MPI_Iprobe (MPI_ANY_SOURCE, ID_MAIL, MPI_COMM_WORLD, &flag, &status);

  if (!flag) return (0);

  *sourceId = status.MPI_SOURCE;
  *reqId    = ID_MAIL;

  if (MPI_SUCCESS != MPI_Recv (data, RAND_MSG_SIZE, MPI_BYTE, *sourceId, ID_MAIL, MPI_COMM_WORLD, &status) )
  {
    printf ("GetMail: %d no mesg!\n", me);
    return (-1);
  }

  MPI_Get_count (&status, MPI_BYTE, size);

  return (1);
}


/****************************************************************************/
/*                                                                          */
/* Miscellaneous                                                            */
/*                                                                          */
/****************************************************************************/

int UsedSpace (void)

{
  return ((int)(100.0*((float)vc_count)/((float)MAXVCHAN)));
}

void PrintHostMessage (char *s)

{
  printf ("%s", s);
}

double CurrentTime (void)

{
  return(((float)(clock())/((float)CLOCKS_PER_SEC)));
}

int Distance (int p, int q)

{
  int pX,pY,pZ,qX,qY,qZ;

  pX = p%DimX;
  p = p/DimX;
  pY = p%DimY;
  pZ = p/DimY;

  qX = q%DimX;
  q = q/DimX;
  qY = q%DimY;
  qZ = q/DimY;

  return (ABS (pX-qX) + ABS (pY-qY) + ABS (pZ-qZ) );
}
