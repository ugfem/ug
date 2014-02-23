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
/*            14 Sep 1995, MPI version                                      */
/*            29 Jan 2003, pV3 concentrator support                         */
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
#include "config.h"
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


#define TRUE  1
#define FALSE 0

/* MPI library */
/*
   #include <mpp/mpi.h>
 */
#include <mpi.h>

#include "../ppif_general.h"
#include "../ppif.h"

/*#include "compiler.h"*/

/*#define _PV3*/

#ifdef _PV3
#include <pV3.h>
#endif

USING_PPIF_NAMESPACE

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
#define PTYPE_ANY       -1L     /* process type: any process                */

#define ID_ARRAY        100     /* channel id: array                        */
#define ID_TREE         101     /* channel id: tree                         */
#define ID_GLOBAL       102     /* channel id: global                       */
#define ID_MAIL         103     /* channel id: mail                         */

#define ABS(i)          (((i)<0) ? (-(i)) : (i))

#define XPOS(aid)       (aid&0xFF)              /* xpos from compact form   */
#define YPOS(aid)       ((aid&0xFF00)>>8)       /* ypos from compact form   */
#define ZPOS(aid)       ((aid&0xFF0000)>>16)    /* zpos from compact form   */

#define PPIF_SUCCESS    0       /* Return value for success                 */
#define PPIF_FAILURE    1       /* Return value for failure                 */

#ifndef _PV3
#define COMM MPI_COMM_WORLD
#else
#define COMM Comm
#endif

/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

typedef struct {
  int p;
  int chanid;
} MPIVChannel;

typedef MPIVChannel *MPIVChannelPtr;


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
int PPIF_NS_PREFIX me;                          /* my processor id                          */
int PPIF_NS_PREFIX master;                      /* id of master processor                   */
int PPIF_NS_PREFIX procs;                       /* number of processors in the network      */

/* 3D array dimensions, may be 1 !          */
int PPIF_NS_PREFIX DimX, PPIF_NS_PREFIX DimY, PPIF_NS_PREFIX DimZ;

/* Tree structure */
int PPIF_NS_PREFIX degree;                      /* degree of downtree nodes                 */
VChannelPtr PPIF_NS_PREFIX uptree;              /* channel uptree                           */
VChannelPtr PPIF_NS_PREFIX downtree[MAXT];      /* channels downtree (may be empty)         */
int PPIF_NS_PREFIX slvcnt[MAXT];                /* number of processors in subtree          */

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static int vc_count=0;          /* number of used VChan                     */

#ifdef _PV3
static MPI_Comm Comm;
#endif

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

int SendSync (void* v, void *data, int size);
int RecvSync (void* v, void *data, int size);

/****************************************************************************/
/*                                                                          */
/* routines for handling virtual channels                                   */
/*                                                                          */
/****************************************************************************/

static VChannelPtr NewVChan (int p, int id)

{
  MPIVChannelPtr myChan = (MPIVChannelPtr)malloc(sizeof(MPIVChannel));

  myChan->p      = p;
  myChan->chanid = id;

  vc_count += 1;

  return (myChan);
}


static void DeleteVChan (VChannelPtr myChan)

{
  free(myChan);

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

int PPIF_NS_PREFIX aid_to_pid (int x, int y, int z)

{
  if ((x<0)||(x>=DimX)) return (-1);
  if ((y<0)||(y>=DimY)) return (-1);
  if ((z<0)||(z>=DimZ)) return (-1);

  return ( (z*DimY+y)*DimX+x);
}

int PPIF_NS_PREFIX pid_to_aid (int p)

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

static void Factor (int N, int *pn, int *pm)

{
  int n = (int)ceil (sqrt ((double) N));
  int m = (int)floor (sqrt ((double) N));

  while (n*m != N)
  {
    if (n*m < N) n++;
    else m--;
  }

  *pn = n; *pm = m;
}



static int PPIFBeganMPI=0; /* remember that PPIF started MPI */


int PPIF_NS_PREFIX InitPPIF (int *argcp, char ***argvp)
{
  int i, succ, sonr, sonl;
  MPI_Status status;
  int mpierror, mpiinitialized;

  /* the following is due to Klaus-Dieter Oertel, 961016;
         (original idea from the developers of the PetSc library)  */
  /* ppif checks whether MPI has been started by another
     library and starts it only if necessary. */
  mpierror = MPI_Initialized(&mpiinitialized);
  if (mpierror) MPI_Abort(MPI_COMM_WORLD, mpierror);
  if (!mpiinitialized)
  {
    mpierror = MPI_Init (argcp, argvp);
    if (mpierror) MPI_Abort( MPI_COMM_WORLD, mpierror);
    PPIFBeganMPI = 1;
  }
#ifdef _PV3
  else
  {
    printf("MPI already initialized, InitPPIF() failed.\n");
    return PPIF_FAILURE;
  }
  if (pV_MPIStart(MPI_COMM_WORLD, 1, 0, 0, &Comm) != 0) {
    printf("pV3 Concentrator cannot be selected. InitPPIF() failed.\n");
    return PPIF_FAILURE;
  }
#endif
  MPI_Comm_rank (COMM, &me);
  MPI_Comm_size (COMM, &procs);

  master = 0;

  vc_count = 0;

  DimZ = 1;
  Factor(procs, &DimX, &DimY);

#ifndef FOR_DUNE
  if (me==master) printf("DimX=%d, DimY=%d, DimZ=%d\n", DimX, DimY, DimZ);
#endif

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
    MPI_Recv ((void *) &(slvcnt[i]), (int) sizeof(int), MPI_BYTE,
              ((MPIVChannel*)downtree[i])->p, ID_TREE, COMM, &status);
    succ += slvcnt[i];
  }
  if (me>0)
  {
    MPI_Send ((void *) &succ, (int) sizeof(int), MPI_BYTE, (int)(me-1)/2, ID_TREE, COMM);
  }

  return (PPIF_SUCCESS);
}


int PPIF_NS_PREFIX ExitPPIF ()
{
  int mpierror;

  if (PPIFBeganMPI)
  {
#ifdef _PV3
    pV_MPISTOP();
#endif
    mpierror = MPI_Finalize();
    if (mpierror) MPI_Abort(MPI_COMM_WORLD, mpierror);
    PPIFBeganMPI = 0;
  }

  return PPIF_SUCCESS;
}


/****************************************************************************/
/*                                                                          */
/* Tree oriented functions                                                  */
/*                                                                          */
/****************************************************************************/

int PPIF_NS_PREFIX Broadcast (void *data, int size)

{
  if (MPI_SUCCESS != MPI_Bcast (data, size, MPI_BYTE, master, COMM) )
    return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int PPIF_NS_PREFIX Concentrate (void *data, int size)

{
  if (me != master)
    if (SendSync (uptree, data, size) < 0) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int PPIF_NS_PREFIX GetConcentrate (int slave, void *data, int size)

{
  if (slave < degree)
    if (RecvSync (downtree[slave], data, size) < 0) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int PPIF_NS_PREFIX Spread (int slave, void *data, int size)

{
  if (slave < degree)
    if (SendSync (downtree[slave], data, size) < 0) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int PPIF_NS_PREFIX GetSpread (void *data, int size)

{
  if (me!=master)
    if (RecvSync (uptree, data, size) < 0) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}

int PPIF_NS_PREFIX Synchronize ()

{
  if (MPI_SUCCESS != MPI_Barrier (COMM) ) return (PPIF_FAILURE);

  return (PPIF_SUCCESS);
}


/****************************************************************************/
/*                                                                          */
/* Synchronous communication                                                */
/*                                                                          */
/****************************************************************************/

VChannelPtr PPIF_NS_PREFIX ConnSync (int p, int id)

{
  return (NewVChan (p, id) );
}

int PPIF_NS_PREFIX DiscSync (void* v)

{
  VChannelPtr vc = (VChannelPtr)v;
  DeleteVChan (vc);

  return (0);
}

int PPIF_NS_PREFIX SendSync (void* v, void *data, int size)

{
  if (MPI_SUCCESS == MPI_Ssend (data, size, MPI_BYTE,
                                ((MPIVChannel*)v)->p, ((MPIVChannel*)v)->chanid, COMM) )
    return (size);
  else
    return (-1);
}

int PPIF_NS_PREFIX RecvSync (void* v, void *data, int size)

{
  int count = -1;
  MPI_Status status;

  if (MPI_SUCCESS == MPI_Recv (data, size, MPI_BYTE,
                               ((MPIVChannel*)v)->p, ((MPIVChannel*)v)->chanid, COMM, &status) )
    MPI_Get_count (&status, MPI_BYTE, &count);

  return (count);
}


/****************************************************************************/
/*                                                                          */
/* Asynchronous communication                                               */
/*                                                                          */
/****************************************************************************/

VChannelPtr PPIF_NS_PREFIX ConnASync (int p, int id)

{
  return (NewVChan (p,id) );
}

int PPIF_NS_PREFIX InfoAConn (void* v)

{
  return (v ? 1 : -1);
}


int PPIF_NS_PREFIX DiscASync (void* v)

{
  DeleteVChan ((MPIVChannel*)v);
  return (PPIF_SUCCESS);
}

int PPIF_NS_PREFIX InfoADisc (void* v)

{
  return (TRUE);
}

#define REQUEST_HEAP

msgid PPIF_NS_PREFIX SendASync (void* v, void *data, int size, int *error)

{
#  ifdef REQUEST_HEAP
  MPI_Request *req;

  if (req = (MPI_Request*)malloc (sizeof (MPI_Request) ) )
  {
    if (MPI_SUCCESS == MPI_Isend (data, size, MPI_BYTE,
                                  ((MPIVChannel*)v)->p, ((MPIVChannel*)v)->chanid, COMM, req) )
    {
      *error = FALSE;
      return ((msgid) req);
    }
  }

#  else
  MPI_Request Req;
  if (MPI_SUCCESS == MPI_Isend (data, size, MPI_BYTE,
                                ((MPIVChannel*)v)->p, ((MPIVChannel*)v)->chanid, COMM, &Req) )
  {
    *error = FALSE;
    return ((msgid) Req);
  }
#  endif

  *error = TRUE;
  return NULL;
}


msgid PPIF_NS_PREFIX RecvASync (void* v, void *data, int size, int *error)

{
#  ifdef REQUEST_HEAP
  MPI_Request *req;

  if (req = (MPI_Request*)malloc (sizeof (MPI_Request) ) )
  {
    if (MPI_SUCCESS == MPI_Irecv (data, size, MPI_BYTE,
                                  ((MPIVChannel*)v)->p, ((MPIVChannel*)v)->chanid, COMM, req) )
    {
      *error = FALSE;
      return ((msgid) req);
    }
  }

#  else
  MPI_Request Req;

  if (MPI_SUCCESS == MPI_Irecv (data, size, MPI_BYTE,
                                ((MPIVChannel*)v)->p, ((MPIVChannel*)v)->chanid, COMM, &Req) )
  {
    *error = FALSE;
    return ((msgid) Req);
  }

#  endif

  *error = TRUE;
  return (NULL);
}


int PPIF_NS_PREFIX InfoASend (void* v, msgid m)

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


int PPIF_NS_PREFIX InfoARecv (void* v, msgid m)
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

int PPIF_NS_PREFIX SendMail (int destId, int reqId, void *data, int size)

{
  if (MPI_SUCCESS == MPI_Send (data, size, MPI_BYTE, destId, ID_MAIL, COMM) )
    return (PPIF_SUCCESS);

  return (PPIF_FAILURE);
}

int PPIF_NS_PREFIX GetMail (int *sourceId, int *reqId, void *data, int *size)

{
  MPI_Status status;
  int flag;

  MPI_Iprobe (MPI_ANY_SOURCE, ID_MAIL, COMM, &flag, &status);

  if (!flag) return (0);

  *sourceId = status.MPI_SOURCE;
  *reqId    = ID_MAIL;

  if (MPI_SUCCESS != MPI_Recv (data, RAND_MSG_SIZE, MPI_BYTE, *sourceId, ID_MAIL, COMM, &status) )
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

int PPIF_NS_PREFIX UsedSpace ()

{
  return vc_count*sizeof(MPIVChannel);
}

void PPIF_NS_PREFIX PrintHostMessage (const char *s)

{
  printf ("%s", s);
}

double PPIF_NS_PREFIX CurrentTime ()

{
  return(((float)(clock())/((float)CLOCKS_PER_SEC)));
}

int PPIF_NS_PREFIX Distance (int p, int q)

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
