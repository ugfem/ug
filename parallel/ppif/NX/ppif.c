// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ppif.c                                                        */
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
/*            Paragon OSF/1 Implementation                                  */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 3a                                                */
/*            7000 Stuttgart 80                                             */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   17 Aug 1992, begin                                            */
/*            18 Feb 1993, Indigo version                                   */
/*            02 Jun 1993, Paragon OSF/1 version                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/* uncomment this to get information about the physical partition shape
   and size on startup (InitPPIF) */
/*#define PPIF_SHOW_PARTITION_INFO*/


/* use the following define to specify how many processes should run
   on one processor. normal setting is 1 process per processor, set
   PROCESS_BITS to 0 for this case. the number of processes per proc is
   2^PROCESS_BITS. */
#define PROCESS_BITS  0


/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

/* standard C library */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>



/* Paragon OSF library */
#include <nx.h>
#include <errno.h>

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

#define RAND_MSG_SIZE 128   /* max size of random messages					*/

#define MAXT        15      /* maximum number of downtree nodes max log2(P) */

#define MAXVCHAN    1024        /* maximum number of channels                                   */

#define PTYPE_ANY   -1L         /* process type: any process					*/

#define ID_ARRAY        100L    /* channel id: array							*/
#define ID_TREE         101L    /* channel id: tree								*/
#define ID_SYNC         102L    /* channel id: synchronize						*/
#define ID_GLOBAL       103L    /* channel id: global							*/
#define ID_MAIL         104L    /* channel id: mail								*/

#define MID_NOMESG      -1L             /* empty message								*/
#define MID_ERROR       -1L             /* error, no message							*/


#define ABS(i) (((i)<0) ? (-(i)) : (i))

#define CHANIDHELP(id,i,j)      (((i)<<20)|((j)<<12)|((id)<<4))
#define CHANID(id,i,j)          ((long)(((i)<(j)) ? CHANIDHELP(id,j,i) : CHANIDHELP(id,i,j)))

#define XPOS(aid)           (aid&0xFF)           /* xpos from compact form  */
#define YPOS(aid)           ((aid&0xFF00)>>8)    /* ypos from compact form  */
#define ZPOS(aid)           ((aid&0xFF0000)>>16) /* zpos from compact form  */




/* macros for mapping PPIF-node numbers to NX-node numbers */
#define PROCESSES_PER_PROC   (1<<PROCESS_BITS)
#define PTYPE_MASK           (PROCESSES_PER_PROC-1)
#define PPIF_PROC(n,t)       (((n)<<PROCESS_BITS)|(t))
#define NX_PROC(p)           ((p)>>PROCESS_BITS)
#if (PROCESS_BITS == 0)
#define NX_PTYPE_ANY(p)      PTYPE_ANY
#define NX_PTYPE(p)          0L
#else
#define NX_PTYPE_ANY(p)      ((p) & PTYPE_MASK)
#define NX_PTYPE(p)          ((p) & PTYPE_MASK)
#endif


/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

enum directions {north,east,south,west,up,down};


typedef struct {
  int used;
  long p;
  long chanid;
  int scnt;                     /* count for preserving msg order */
  int rcnt;                     /* count for preserving msg order */
} VChannel;

typedef VChannel *VChannelPtr;
typedef int msgid;


/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
RCSID("$Header$",PPIF_RCS_STRING)


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

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
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static VChannel *VChan;
static int vc_count=0;            /* number of used VChan                   */
static int vc_free=0;             /* rotating pointer to find free VChan    */

static char buf[100];    /* string buffer for debugging */


/****************************************************************************/
/*                                                                          */
/* routines for handling virtual channels                                                           */
/*                                                                          */
/****************************************************************************/

int InitVChan (void)
{
  int i;

  VChan = (VChannel *) malloc(sizeof(VChannel)*MAXVCHAN);
  if (VChan==0L)
    return(0);

  for (i=0; i<MAXVCHAN; i++) VChan[i].used = FALSE;

  vc_count=0;
  vc_free=0;

  return(1);
}


VChannelPtr NewVChan (int p, int id)
{
  int i;
  VChannelPtr myChan;
  if (vc_count<MAXVCHAN)
  {
    while (VChan[vc_free].used) {
      if (++vc_free==MAXVCHAN)
        vc_free=0;
    }
    myChan = &VChan[vc_free];
    myChan->used = TRUE;
    myChan->p  = p;
    myChan->chanid = CHANID((long)id,(long)p,(long)me);
    myChan->rcnt = 0;
    myChan->scnt = 0;

    vc_count++;
    vc_free = (vc_free+1)%MAXVCHAN;
    return(myChan);
  }
  else
  {
    printf("%4d: PPIF-error, no more VChannels in NewVChan(), "
           "dest=%d, id=%d\n", me, p, id);
    return(NULL);
  }
}


void DeleteVChan (VChannelPtr myChan)
{
  myChan->used = FALSE;
  vc_count--;
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

void Factor(N, pn, pm)
int N, *pn, *pm;
/*
   Factor N into two integers that are as close together as possible
 */
{
  int n = ceil(sqrt((double) N));
  int m = floor(sqrt((double) N));

  while (n*m != N) {
    if (n*m < N)
      n++;
    else
      m--;
  }
  *pn = n; *pm = m;
}


int InitPPIF (int *argcp, char ***argvp)
{
  int forks, i, succ, sonr, sonl, aid;
  long local_node = mynode();
  long pid_list[PROCESSES_PER_PROC];
  int child_node_offset = 0;

  for(forks=1; child_node_offset==0 && forks<PROCESSES_PER_PROC; forks++)
  {
    long child_spawned =
      nx_nfork(&local_node, 1L, (long)forks, &(pid_list[forks]));

    if (child_spawned==-1)
    {
      fprintf(stderr, "%4d: PPIF: couldn't spawn child process (errno=%ld)\n",
              local_node, (long)errno);
      exit(1);
    }

    if (child_spawned==0)
    {
      /* I am a successfully spawned child! */
      child_node_offset = forks;
    }
  }

  me = PPIF_PROC(mynode(), child_node_offset);
  master = 0;
  procs = numnodes()<<PROCESS_BITS;


  if (myptype()!=(long)child_node_offset) {
    setptype((long)child_node_offset);
    printf("%4d: PPIF-warning: corrected PTYPE to %d\n",
           me, child_node_offset);
  }

  if (!InitVChan())
  {
    printf("%4d: PPIF-error: couldn't get VChannel memory!\n", me);
  }

#ifdef PPIF_SHOW_PARTITION_INFO
  if (me==0) {
    nx_nodes_t mynodes;
    unsigned long nnodes;
    int i, status;
    long rows, cols;

    status = nx_app_nodes(0, &mynodes, &nnodes);

    assert(nnodes==procs || PROCESSES_PER_PROC>1);

    if (status!=0)
    {
      printf("%4d: PPIF-warning: couldn't get partition-info\n", me);
    }
    else
    {
      printf("PARTITION_INFO:\n");
      for(i = 0; i < nnodes; i++)
      {
        if (i%8==0) printf("\t%4d: ", i);
        printf(" %4d", mynodes[i]);
        if (i%8==7) printf("\n");
      }
      printf("\n");
    }

    status = nx_app_rect(&rows, &cols);
    if (status!=0)
    {
      printf("%4d: PPIF-warning: couldn't get partition-rect\n", me);
    }
    else
    {
      printf("PARTITION_RECT: row=%ld x col=%ld\n", rows, cols);
    }
  }
#endif

  /* 3D array configuration */
  /*
          if (me==master) {
                  do {
                          printf("DimX? "); scanf("%d", &DimX);
                          printf("DimY? "); scanf("%d", &DimY);
                          printf("DimZ? "); scanf("%d", &DimZ);
                          if (DimX*DimY*DimZ != procs)
                                  printf("Procs=%d, nochmal!\n", procs);
                  } while (DimX*DimY*DimZ != procs);
          }
   */
  DimZ = 1;
  Factor(procs, &DimX, &DimY);
  if (me==master)
  {
    fprintf(stderr, "PPIF: DimX=%d, DimY=%d, DimZ=%d\n", DimX, DimY, DimZ);
    if (PROCESSES_PER_PROC!=1)
      fprintf(stderr, "PPIF: ProcessesPerNode=%d\n", PROCESSES_PER_PROC);

  }

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

  if (sonl<procs) {
    degree++;
    downtree[0] = NewVChan(sonl,ID_TREE);
  } else {
    downtree[0] = NULL;
  }
  if (sonr<procs) {
    degree++;
    downtree[1] = NewVChan(sonr,ID_TREE);
  } else {
    downtree[1] = NULL;
  }

  if (me>0) {
    uptree = NewVChan((me-1)/2,ID_TREE);
  } else {
    uptree = NULL;
  }

  succ=1;
  for(i=0; i<degree; i++)
  {
    long p = downtree[i]->p;

    crecvx(ID_GLOBAL, (char *)&(slvcnt[i]), (long)sizeof(int),
           NX_PROC(p), NX_PTYPE_ANY(p), msginfo);
    succ += slvcnt[i];
  }
  if (me>0)
  {
    long p = (long)(me-1)/2;
    csend(ID_GLOBAL, (char *)&succ, (long)sizeof(int),
          NX_PROC(p), NX_PTYPE(p));
  }
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
  long ret;

  if (me==master)
  {
    long ptypes;
    for(ptypes=0; ptypes<PROCESSES_PER_PROC; ptypes++)
    {
      /*printf(" SEND  %3d -> all  mit %lx\n", me, ID_GLOBAL);*/
      ret = _csend(ID_GLOBAL, (char *)data, (long)size, -1L, ptypes);
      if (ret!=0L)
        return(1);
    }
    return(0);
  }
  else
  {
    /*printf(" RECV  sby -> %3d  mit %lx\n", me, ID_GLOBAL);*/
    ret = _crecv(ID_GLOBAL, (char *)data, (long)size);
    return((ret!=-1L) ? 0 : 1);
  }
}

int Concentrate (void *data, int size)
{
  long ret;

  if (me!=master)
  {
    long p = uptree->p;
    ret = _csend(uptree->chanid|uptree->scnt, (char *)data, (long)size,
                 NX_PROC(p), NX_PTYPE(p));
    uptree->scnt = (uptree->scnt+1)&0x0f;
    return((ret==0L) ? 0 : 1);
  }
  return(0);
}

int GetConcentrate (int slave, void *data, int size)
{
  VChannel *chan;
  long ret;

  if (slave<degree)
  {
    long p;
    chan = downtree[slave];
    p = chan->p;
    ret = _crecvx(chan->chanid|chan->rcnt, (char *)data, (long)size,
                  NX_PROC(p), NX_PTYPE_ANY(p), msginfo);
    chan->rcnt = (chan->rcnt+1)&0x0f;
    return((ret==0L) ? 0 : 1);
  }
  return(0);
}

int Spread (int slave, void *data, int size)
{
  VChannel *chan;
  long ret;

  if (slave<degree)
  {
    long p;
    chan = downtree[slave];
    p = chan->p;
    ret = _csend(chan->chanid|chan->scnt, (char *)data, (long)size,
                 NX_PROC(p), NX_PTYPE(p));
    chan->scnt = (chan->scnt+1)&0x0f;
    return((ret==0L) ? 0 : 1);
  }
  return(0);
}

int GetSpread (void *data, int size)
{
  long ret;

  if (me!=master)
  {
    long p = uptree->p;
    ret = _crecvx(uptree->chanid|uptree->rcnt, (char *)data, (long)size,
                  NX_PROC(p), NX_PTYPE_ANY(p), msginfo);
    uptree->rcnt = (uptree->rcnt+1)&0x0f;
    return((ret==0L) ? 0 : 1);
  }
  return(0);
}



int Synchronize (void)
{
#if (PROCESS_BITS == 0)
  long ret;

  ret = _gsync();
  return((ret!=-1L) ? 0 : 1);
#else
  int l, dummy=0;

  for(l=degree-1; l>=0; l--) {
    GetConcentrate(l, &dummy, sizeof(dummy));
  }
  Concentrate(&dummy, sizeof(dummy));

  Broadcast(&dummy, sizeof(dummy));
#endif
}




/****************************************************************************/
/*                                                                          */
/* Synchronous communication                                                */
/*                                                                          */
/****************************************************************************/

VChannelPtr ConnSync (int p, int id)
{
  return(NewVChan(p,id));
}

int DiscSync (VChannelPtr vc)
{
  DeleteVChan(vc);
  return(0);
}

int SendSync (VChannelPtr vc, void *data, int size)
{
  long ret, p = vc->p;

  ret = _csend(vc->chanid, (char *)data, (long)size,
               NX_PROC(p), NX_PTYPE(p));
  return((ret==0) ? size : -1);
}

int RecvSync (VChannelPtr vc, void *data, int size)
{
  long ret, p = vc->p;

  ret = _crecvx(vc->chanid, (char *)data, (long)size,
                NX_PROC(p), NX_PTYPE_ANY(p), msginfo);
  return((int)ret);
}


/****************************************************************************/
/*                                                                          */
/* Asynchronous communication                                               */
/*                                                                          */
/****************************************************************************/

VChannelPtr ConnASync (int p, int id)
{
  return(NewVChan(p,id));
}

int InfoAConn (VChannelPtr vc)
{
  return((vc!=NULL) ? 1 : -1);
}


int DiscASync (VChannelPtr vc)
{
  DeleteVChan(vc);
  return(0);
}

int InfoADisc (VChannelPtr vc)
{
  return(1);
}


msgid SendASync (VChannelPtr vc, void *data, int size, int *error)
{
  long id, p = vc->p;

  *error = 0;
  id = _isend(vc->chanid|vc->scnt, (char *)data, (long)size,
              NX_PROC(p), NX_PTYPE(p));

  /*sprintf(buf, "    %d: SA   id-e=%ld   an %d  vc %8p  id-i %8p  size %4d", me, id, vc->p, vc,vc->chanid|vc->scnt, size);*/
  vc->scnt = (vc->scnt+1)&0x0f;

  if (id==MID_ERROR) {
    /*printf("%s\n", buf);*/
    *error = 1;
    return(0);
  }

  return(id);
}


msgid RecvASync (VChannelPtr vc, void *data, int size, int *error)
{
  long id, p = vc->p;

  *error = 0;
  id = _irecvx(vc->chanid|vc->rcnt, (char *)data, (long)size,
               NX_PROC(p), NX_PTYPE_ANY(p), msginfo);

  /*sprintf(buf, "    %d: RA   id-e=%ld  von %d  vc %8p  id-i %8p  size %4d", me, id, vc->p, vc,vc->chanid|vc->rcnt,size);*/
  vc->rcnt = (vc->rcnt+1)&0x0f;

  if (id==MID_ERROR) {
    /*printf("%s\n", buf);*/
    *error = 1;
    return(0);
  }

  return(id);
}


int InfoASend (VChannelPtr vc, msgid m)
{
  long ret;

  ret = _msgdone((long) m);
  /*printf("    %d: IS   ret=%ld   an %d  vc %8p  id %8p\n", me, ret, vc->p, vc,vc->chanid);*/
  switch (ret) {
  case 1 :                              /* message transfer has been completed		*/
    return(1);

  case 0 :                              /* transfer not yet complete						*/
    return(0);

  default :                             /* error during transfer							*/
    fprintf(stdout,
            "%4d: PPIF-error in InfoASend(), errno=%ld (see errno.h)\n",
            me, (long)errno);
    return(-1);
  }
}


int InfoARecv (VChannelPtr vc, msgid m)
{
  long ret;

  ret = _msgdone((long) m);
  /*printf("    %d: IR   ret=%ld  von %d  vc %8p  id %8p\n", me, ret, vc->p, vc,vc->chanid);*/
  switch (ret) {
  case 1 :                              /* message transfer has been completed		*/
    return(1);

  case 0 :                              /* transfer not yet complete						*/
    return(0);

  default :                             /* error during transfer							*/
    fprintf(stdout,
            "%4d: PPIF-error in InfoARecv(), errno=%ld (see errno.h)\n",
            me, (long)errno);
    return(-1);
  }
  return(0);
}


/****************************************************************************/
/*                                                                          */
/* Random communication                                                         */
/*                                                                          */
/****************************************************************************/

/*
   extern long errno;
 */

int SendMail (int destId, int reqId, void *data, int size)
{
  long ret, p = (long)destId;

  ret = _csend(ID_MAIL, (char *)data, (long)size, NX_PROC(p), NX_PTYPE(p));
  return((ret==0) ? 0 : 1);
}

int GetMail (int *sourceId, int *reqId, void *data, int *size)
{
  long ret;

  ret = _iprobe(ID_MAIL);
  if (ret!=1L) return((int) ret);


  ret = _crecv(ID_MAIL, (char *)data, (long)RAND_MSG_SIZE);
  if (ret==-1L) {
    printf("GetMail: %d no mesg! errno=%d (see errno.h)\n", me, (int)errno);
    return(-1);
  }

  *sourceId = (int)(PPIF_PROC(infonode(), infoptype()));
  if (infonode()==-1)
    *sourceId = -1;

  *reqId = (int)(ID_MAIL);
  *size = (int)infocount();

  return((*sourceId!=-1L)&&(*reqId!=-1L) ? 1 : -1);
}


/****************************************************************************/
/*                                                                          */
/* Miscellaneous                                                                        */
/*                                                                          */
/****************************************************************************/

int UsedSpace (void)
{
  return((int)(100.0*((float)vc_count)/((float)MAXVCHAN)));
}

void PrintHostMessage (char *s)
{
  printf("%s",s);
}

double CurrentTime (void)
{
  return(dclock());
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

  return(ABS(pX-qX)+ABS(pY-qY)+ABS(pZ-qZ));
}
