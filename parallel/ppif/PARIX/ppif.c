// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      parix-ppif.c                                                  */
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
/*            Parsytec PARIX 1.0 Implementation                             */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*            internet: bastian@iwr1.iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   17 August 1992, begin                                         */
/*            95/10/04 kb  added Spread/GetSpread                           */
/*            95/11/06 kb  changed parameters for InitPPIF()                */
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
#include <string.h>

/* PARIX library */
#include <sys/root.h>
#include <sys/sys_rpc.h>
#include <sys/topology.h>
#include <sys/rrouter.h>
#include <sys/logerror.h>
#include <sys/thread.h>
#include <sys/time.h>

/* LocalLink has no definition */
int LocalLink (LinkCB_t *Link[2]);

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define MAXSLOTS    128     /* maximum number of channels in async mode     */
#define ASYNCMAX        512             /* max number of asynchronous messages in fifos */
#define MAILBOXSIZE 128         /* number of random messages per processor		*/

#define MAXT        15      /* maximum number of downtree nodes max log2(P) */
#define ARRAYID     10000   /* channel id's for 3D array structure          */
#define TREEID      10001   /* channel id's for tree structure              */
#define STACKSIZE   8192     /* stack size for threads in async comm         */
#define RMSTACKSIZE 8192     /* stack size random comm receiver process         */

#define RAND_MSG_SIZE 64    /* short random messages                                            */
#define MSGTYPE     1       /* message type number used for random comm		*/

#define VERBOSE         0               /* verbose level								*/

#define DELAY       100     /* deschedule for a number of low prio ticks    */

/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

typedef LinkCB_t *VChannelPtr;
typedef struct messagecarrier *msgid ;

enum directions {north,east,south,west,up,down};

struct fifo {
  struct messagecarrier *first,*last;        /* list of messages */
  Semaphore_t sema;                          /* controls access to the list */
  int count;                                 /* no of messages in fifo */
  VChannelPtr channel[2];                            /* wait if fifo empty */
} ;

typedef struct fifo FIFO;
typedef struct messagecarrier MessageCarrier;

typedef struct {
  int used;                                             /* true if record is used					*/
  VChannelPtr Destination;              /* synchronous channel to destination proc  */
  FIFO outgoing;                                /* holds outgoing message carriers			*/
  FIFO incoming;                                /* holds incoming message carriers			*/
  int SendFlag;                                 /* counts outgoing messages (only 0,1)      */
  int RecvFlag;                                 /* counts incoming messages (only 0,1)      */
  Thread_t *SendThread;             /* outgoing buffer process                  */
  Thread_t *RecvThread;                 /* incoming buffer process                  */
} Slot_t ;

typedef int (*FunctionPtr)(void) ;

typedef struct {
  int used;
  int sourceId;
  int reqId;
  int size;
  char data[RAND_MSG_SIZE];
} RandomMessage ;

struct messagecarrier {
  struct messagecarrier *next;
  short used;
  short cmd;
  int ready;
  int size;
  void *data;
} ; /* 16 Bytes */



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

static Slot_t Slots[MAXSLOTS];    /* control structure for async comm       */
static Semaphore_t mySema;        /* semaphore for Slots array              */
static int usedCount=0;         /* number of used slots                   */
static int slotPtr=0;                     /* rotating pointer to find free slot		*/

static MessageCarrier mc[ASYNCMAX]; /* stores message parameters				*/
static Semaphore_t mcSema;                /* control array access					*/
static int mcPtr=0;                               /* rotating pointer to find free mc		*/
static int mcCount=0;                     /* number of used mc's					*/

static RandomMessage rm[MAILBOXSIZE]; /* store random messages				*/
static Semaphore_t rmSema;                /* access to mail box						*/
static int rmCount=0;                     /* number of messages in box				*/
static int rmRecvPtr=0;                   /* rotating pointer for receiver			*/
static int rmReadPtr=0;                   /* rotating pointer for read				*/
static int rmError=0;                     /* number of messages thrown away			*/
static int rmExit=0;                      /* exit random communication				*/
static Thread_t *rmThread=NULL;   /* random communication process			*/

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

#define ABS(i) (((i)<0) ? (-(i)) : (i))
int GetMailProcess (void);
int SendMail (int destId, int reqId, void *data, int size);

/****************************************************************************/
/*                                                                          */
/* message carrier fifos for asynchronous communication                         */
/*                                                                          */
/****************************************************************************/

static int Fifo_init (FIFO *f)
{
  int error;

  f->first = NULL;
  f->last = NULL;
  f->count = 0;
  InitSem(&(f->sema),1);
  error = LocalLink(&(f->channel[0]));
  if (error<0) return(1);
  return(0);
}


static int Fifo_exit (FIFO *f)
{
  BreakLink(f->channel[0]);
  BreakLink(f->channel[1]);
  return(0);
}


static int Fifo_count (FIFO *f)
{
  return(f->count);
}

static void Fifo_in (FIFO *f,MessageCarrier *item)
{
  Wait(&(f->sema));                      /* wait for access */
  if (f->count<0)                        /* is somebody waiting ? */
  {
    (f->count)++;                            /* message has been sent */
    SendLink(f->channel[0],(void *)(&item),sizeof(item));
  }
  else                                   /* put message in queue */
  {
    if (f->first==NULL)                      /* the FIFO is empty */
      f->first = item;                           /* place first item in list */
    else
      f->last->next = item;                      /* or append item */
    f->last = item;
    item->next = NULL;                       /* and without successor */
    (f->count)++;                            /* message has been sent */
  }
  Signal(&(f->sema));                    /* leave critical section */
}



static MessageCarrier *Fifo_out (FIFO *f)
{
  MessageCarrier *item;

  Wait(&(f->sema));                      /* enter critical section */
  if (f->count==0)                       /* queue is empty */
  {
    (f->count)--;                            /* I am waiting */
    Signal(&(f->sema));                      /* allow others to come in */
    RecvLink(f->channel[1],&item,sizeof(item));
  }
  else
  {
    item = f->first;                         /* remove first item in list */
    f->first = item->next;                   /* update first */
    (f->count)--;
    Signal(&(f->sema));
  }
  return(item);
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


static int BuildTree (int start, int len, int dim)
{
  int error,mypos,l1,l2;

  /* len = 1, do nothing */
  if (len==1) return(0);

  /* get my position */
  switch (dim)
  {
  case 1 : mypos = MyX; break;
  case 2 : mypos = MyY; break;
  case 3 : mypos = MyZ; break;
  default : return(1);
  }

  /* downtree connections */
  if (((len==2)&&(mypos==start))||((len==3)&&(mypos<=start+1)))
  {
    switch (dim)
    {
    case 1 : downtree[degree] = ConnectLink(aid_to_pid(MyX+1,MyY,MyZ),TREEID,&error); break;
    case 2 : downtree[degree] = ConnectLink(aid_to_pid(MyX,MyY+1,MyZ),TREEID,&error); break;
    case 3 : downtree[degree] = ConnectLink(aid_to_pid(MyX,MyY,MyZ+1),TREEID,&error); break;
    default : return(1);
    }
    if (downtree[degree]==NULL) return(1);
    degree++;
  }

  /* uptree connections */
  if (((len==2)&&(mypos==start+1))||((len==3)&&(mypos>=start+1)))
  {
    switch (dim)
    {
    case 1 : uptree = ConnectLink(aid_to_pid(MyX-1,MyY,MyZ),TREEID,&error); break;
    case 2 : uptree = ConnectLink(aid_to_pid(MyX,MyY-1,MyZ),TREEID,&error); break;
    case 3 : uptree = ConnectLink(aid_to_pid(MyX,MyY,MyZ-1),TREEID,&error); break;
    default : return(1);
    }
    if (uptree==NULL) return(1);
  }

  /* recursive subdivision */
  if (len>3)
  {
    l1 = len/2;
    l2 = len-l1;
    if (mypos<start+l1)
    {
      if (mypos==start)
      {
        switch (dim)
        {
        case 1 : downtree[degree] = ConnectLink(aid_to_pid(start+l1,MyY,MyZ),TREEID,&error); break;
        case 2 : downtree[degree] = ConnectLink(aid_to_pid(MyX,start+l1,MyZ),TREEID,&error); break;
        case 3 : downtree[degree] = ConnectLink(aid_to_pid(MyX,MyY,start+l1),TREEID,&error); break;
        default : return(1);
        }
        if (downtree[degree]==NULL) return(1);
        degree++;
      }
      return(BuildTree(start,l1,dim));
    }
    else
    {
      if (mypos==start+l1)
      {
        switch (dim)
        {
        case 1 : uptree = ConnectLink(aid_to_pid(start,MyY,MyZ),TREEID,&error); break;
        case 2 : uptree = ConnectLink(aid_to_pid(MyX,start,MyZ),TREEID,&error); break;
        case 3 : uptree = ConnectLink(aid_to_pid(MyX,MyY,start),TREEID,&error); break;
        default : return(1);
        }
        if (uptree==NULL) return(1);
      }
      return(BuildTree(start+l1,l2,dim));
    }
  }

  /* return ok */
  return(0);
}


int InitPPIF (int *argcp, char ***argvp)
{
  RootProc_t *r;
  int i,j,error,cnt;

  /* get pointer to root structure */
  r = GET_ROOT()->ProcRoot;

  /* set id's */
  me = (r->MyZ*r->DimY + r->MyY)*r->DimX + r->MyX;
  master = 0;
  procs = r->nProcs;

  /* 3D array configuration */
  MyX = r->MyX;
  MyY = r->MyY;
  MyZ = r->MyZ;
  DimX = r->DimX;
  DimY = r->DimY;
  DimZ = r->DimZ;
  for (i=0; i<6; i++) nn[i] = NULL;
  /* alloc east-west */
  if (DimX>1)
  {
    if (MyX>0)
    {
      nn[west] = ConnectLink(me-1,ARRAYID,&error);
      if (nn[west]==NULL) return(1);
    }
    if (MyX<DimX-1)
    {
      nn[east] = ConnectLink(me+1,ARRAYID,&error);
      if (nn[east]==NULL) return(1);
    }
  }
  /* alloc north-south */
  if (DimY>1)
  {
    if (MyY>0)
    {
      nn[south] = ConnectLink(me-DimX,ARRAYID,&error);
      if (nn[south]==NULL) return(1);
    }
    if (MyY<DimY-1)
    {
      nn[north] = ConnectLink(me+DimX,ARRAYID,&error);
      if (nn[north]==NULL) return(1);
    }
  }
  /* alloc up-down */
  if (DimZ>1)
  {
    if (MyZ>0)
    {
      nn[down] = ConnectLink(me-DimX*DimY,ARRAYID,&error);
      if (nn[down]==NULL) return(1);
    }
    if (MyZ<DimZ-1)
    {
      nn[up] = ConnectLink(me+DimX*DimY,ARRAYID,&error);
      if (nn[up]==NULL) return(1);
    }
  }

  /* tree configuration */
  degree = 0;
  uptree = NULL;
  for (i=0; i<MAXT; i++) downtree[i] = NULL;
  /* allocate uptree channel */
  if ((MyX==0)&&(MyY==0))
  {
    error = BuildTree(0,DimZ,3);
    if (error!=0)
    {
      printe("%d: error in BuildTree z-direction\n",me);
      return(1);
    }
  }
  if (MyX==0)
  {
    error = BuildTree(0,DimY,2);
    if (error!=0)
    {
      printe("%d: error in BuildTree y-direction\n",me);
      return(1);
    }
  }
  error = BuildTree(0,DimX,1);
  if (error!=0)
  {
    printe("%d: error in BuildTree x-direction\n",me);
    return(1);
  }

  /* count subtree nodes */
  for (i=0; i<MAXT; i++) slvcnt[i] = 0;
  for (i=0; i<degree; i++)
  {
    if (RecvLink(downtree[i],&j,sizeof(int))!=sizeof(int)) return(1);
    slvcnt[i] = j;
  }
  cnt = 1;
  for (i=0; i<degree; i++) cnt += slvcnt[i];
  if (me!=master)
    if (SendLink(uptree,&cnt,sizeof(int))!=sizeof(int)) return(1);

  /* init structures for arbitrary communication network */
  InitSem(&mySema,1);
  InitSem(&mcSema,1);
  for (i=0; i<MAXSLOTS; i++) Slots[i].used = 0;
  for (i=0; i<ASYNCMAX; i++) mc[i].used = 0;
  for (i=0; i<ASYNCMAX; i++) mc[i].cmd = 0;

  /* initialize random communication */
  InitSem(&rmSema,1);
  for (i=0; i<MAILBOXSIZE; i++) rm[i].used = 0;
  if ((rmThread=StartThread(GetMailProcess,RMSTACKSIZE,&error,0))==NULL) return(1);

  return(0);
}

void ExitPPIF (void)
{
  int dummy;

  rmExit = 1;
  SendMail(me,1,&dummy,sizeof(int));
  WaitThread(rmThread,&dummy);

  return;
}

/****************************************************************************/
/*                                                                          */
/* Tree oriented functions                                                  */
/*                                                                          */
/****************************************************************************/

int Broadcast (void *data, int size)
{
  int i,error;

  error = 0;

  if (me!=master)
    if (RecvLink(uptree,data,size)!=size) error = 1;

  for (i=0; i<degree; i++)
    if (SendLink(downtree[i],data,size)!=size) error = 1;

  return(error);
}

int Concentrate (void *data, int size)
{
  if (me!=master)
    if (SendLink(uptree,data,size)!=size) return(1);
  return(0);
}

int GetConcentrate (int slave, void *data, int size)
{
  if (slave<degree)
    if (RecvLink(downtree[slave],data,size)!=size) return(1);
  return(0);
}


int Spread (int slave, void *data, int size)
{
  if (slave<degree)
    if (SendLink(downtree[slave],data,size)!=size) return(1);
  return(0);
}

int GetSpread (void *data, int size)
{
  if (me!=master)
    if (RecvLink(uptree,data,size)!=size) return(1);
  return(0);
}


int Synchronize (void)
{
  int i,n;

  for (i=degree-1; i>=0; i--) RecvLink(downtree[i],&n,sizeof(int));

  if (me!=master)
  {
    SendLink(uptree,&n,sizeof(int));
    RecvLink(uptree,&n,sizeof(int));
  }

  for (i=0; i<degree; i++) SendLink(downtree[i],&n,sizeof(int));

  return(0);
}


/****************************************************************************/
/*                                                                          */
/* Synchronous communication                                                */
/*                                                                          */
/****************************************************************************/

VChannelPtr ConnSync (int p, int id)
{
  int error;

  return(ConnectLink(p,id,&error));
}

int DiscSync (VChannelPtr vc)
{
  return(-BreakLink(vc));
}

int SendSync (VChannelPtr vc, void *data, int size)
{
  return(SendLink(vc,data,size));
}

int RecvSync (VChannelPtr vc, void *data, int size)
{
  return(RecvLink(vc,data,size));
}


/****************************************************************************/
/*                                                                          */
/* Asynchronous communication                                               */
/*                                                                          */
/****************************************************************************/

static int SenderThread (int p, int id, Slot_t *mySlot)
{
  int error;
  VChannelPtr toDestination;
  MessageCarrier *mess;

  /* make me a high priority process */
  ChangePriority(HIGH_PRIORITY);

  /* connect to destination (sender only) */
  toDestination = ConnectLink(p,id,&error);
  if (toDestination==NULL)
  {
    printe("%d: error in SenderThread() during ConnectLink(), dest=%d, id=%d, error=%d\n",me,p,id,error);
    mySlot->SendFlag = -2;
    return(-2);
  }
  mySlot->Destination = toDestination;

  /* set ready flag */
  mySlot->SendFlag = 0;

  if (VERBOSE>=1)
  {
    printe("%d: SThread started vc=%x p=%d id=%d slot=%x sf=%d rf=%d\n",me,toDestination,p,id,mySlot,mySlot->SendFlag,mySlot->RecvFlag);
  }

  /* send loop */
  while (1)
  {
    /* get next message from fifo */
    mess = Fifo_out(&(mySlot->outgoing));
    if (mess->cmd)
    {
      mySlot->SendFlag = 1;
      break;                   /* break link */
    }

    if (VERBOSE>=1)
    {
      printe("%d: SThread bef send vc=%x p=%d id=%d s=%d data=%x fifo=%d\n",me,toDestination,p,id,mess->size,mess->data,Fifo_count(&(mySlot->outgoing)));
    }

    /* send to destination */
    if ((error=SendLink(toDestination,mess->data,mess->size))<0)
    {
      printe("%d: error in SenderThread() during SendLink(), dest=%d, id=%d size=%d\n",me,p,id,mess->size);
      mySlot->SendFlag = -4;
    }
    /* set ready flag */
    mess->ready = 1;

    if (VERBOSE>=1)
    {
      printe("%d: SThread aft send vc=%x p=%d id=%d s=%d data=%x fifo=%d\n",me,toDestination,p,id,mess->size,mess->data,Fifo_count(&(mySlot->outgoing)));
    }
  }

  /* make me a low priority process */
  ChangePriority(LOW_PRIORITY);

  /* disconnect from destination (sender only) */
  error = BreakLink(toDestination);
  if (error<0)
  {
    printe("%d: error in SenderThread() during 1. BreakLink(), dest=%d, id=%d, error=%d\n",me,p,id,error);
    return(1);
  }

  if (VERBOSE>=2)
  {
    printe("%d: SThread EXITING vc=%x p=%d id=%d slot=%x\n",me,toDestination,p,id,mySlot);
  }

  /* return ok */
  return(0);
}


static int ReceiverThread (int p, int id, Slot_t *mySlot)
{
  int error;
  MessageCarrier *mess;
  VChannelPtr *fromDestination;

  /* make me a high priority process */
  ChangePriority(HIGH_PRIORITY);

  /* set ready flag */
  fromDestination = &(mySlot->Destination);
  mySlot->RecvFlag = 0;

  if (VERBOSE>=1)
  {
    printe("%d: RThread started vc=%x p=%d id=%d slot=%x sf=%d rf=%d\n" ,me,*fromDestination,p,id,mySlot,mySlot->SendFlag,mySlot->RecvFlag);
  }

  /* receive loop */
  while (1)
  {
    /* get next message from fifo */
    mess = Fifo_out(&(mySlot->incoming));
    if (mess->cmd)
    {
      mySlot->RecvFlag = 1;
      break;                   /* break link */
    }

    if (VERBOSE>=1)
    {
      printe("%d: RThread bef recv vc=%x p=%d id=%d s=%d data=%x fifo=%d\n",me,*fromDestination,p,id,mess->size,mess->data,Fifo_count(&(mySlot->incoming)));
    }

    /* get from destination */
    if ((error=RecvLink(*fromDestination,mess->data,mess->size))<0)
    {
      printe("%d: error in ReceiverThread() during RecvLink(), src=%d, id=%d error=%d size=%d\n",me,p,id,error,mess->size);
      mySlot->RecvFlag = -4;
    }
    if (VERBOSE>=1)
    {
      printe("%d: RThread aft recv vc=%x p=%d id=%d s=%d data=%x fifo=%d\n",me,*fromDestination,p,id,mess->size,mess->data,Fifo_count(&(mySlot->incoming)));
    }
    /* set ready flag */
    mess->ready = 1;
  }

  /* make me a low priority process */
  ChangePriority(LOW_PRIORITY);

  if (VERBOSE>=2)
  {
    printe("%d: RThread EXITING vc=%x p=%d id=%d slot=%x\n" ,me,*fromDestination,p,id,mySlot);
  }

  /* return ok */
  return(0);
}

msgid SendASync (VChannelPtr vc, void *data, int size, int *error)
{
  MessageCarrier *mess;
  Slot_t *mySlot;
  int i;

  *error = 0;

  /* cast to slot pointer */
  mySlot = (Slot_t *) vc;

  /* allocate a free message carrier */
  for (i=0; i<ASYNCMAX; i++)
    if (mc[mcPtr].used)
      mcPtr = (mcPtr+1)%ASYNCMAX;
    else
      break;

  /* if no free carrier return error */
  if (i>=ASYNCMAX) {*error = 1; return(NULL);}
  mess = &(mc[mcPtr]);
  mess->used = 1;
  mcPtr = (mcPtr+1)%ASYNCMAX;

  /* put message in fifo */
  mess->ready = 0;
  mess->data = data;
  mess->size = size;
  Fifo_in(&(mySlot->outgoing),mess);
  if (VERBOSE>=2)
  {
    printe("%d: SendASync put a msg, vc=%x, carrier=%x, size=%d, data=%x, fifo=%d\n",me,mySlot,mess,mess->size,mess->data,Fifo_count(&(mySlot->outgoing)));
  }

  /* ok */
  return(mess);
}

msgid RecvASync (VChannelPtr vc, void *data, int size, int *error)
{
  MessageCarrier *mess;
  Slot_t *mySlot;
  int i;

  *error = 0;

  /* cast to slot pointer */
  mySlot = (Slot_t *) vc;

  /* allocate a free message carrier */
  for (i=0; i<ASYNCMAX; i++)
    if (mc[mcPtr].used)
      mcPtr = (mcPtr+1)%ASYNCMAX;
    else
      break;
  if (i>=ASYNCMAX) {*error = 1; return(NULL);}
  mess = &(mc[mcPtr]);
  mess->used = 1;
  mcPtr = (mcPtr+1)%ASYNCMAX;

  /* put message in fifo */
  mess->ready = 0;
  mess->data = data;
  mess->size = size;
  Fifo_in(&(mySlot->incoming),mess);
  if (VERBOSE>=2)
  {
    printe("%d: RecvASync put a msg, vc=%x, carrier=%x, size=%d, data=%x, fifo=%d\n",me,mySlot,mess,mess->size,mess->data,Fifo_count(&(mySlot->incoming)));
  }

  return(mess);
}

int InfoASend (VChannelPtr vc, msgid m)
{
  Slot_t *mySlot;
  int i;

  /* cast to slot pointer */
  mySlot = (Slot_t *) vc;

  /* error in communication */
  if (mySlot->SendFlag<0) return(-1);
  if (m->used==0) return(-1);

  /* if message not ready, deschedule and return 0 */
  if (!m->ready)
  {
    TimeWait(TimeNow()+DELAY);
    return(0);
  }

  /* message arrived, recycle message carrier and return */
  m->used = 0;
  return(1);
}

int InfoARecv (VChannelPtr vc, msgid m)
{
  Slot_t *mySlot;
  int i;

  /* cast to slot pointer */
  mySlot = (Slot_t *) vc;

  /* check error */
  if (mySlot->RecvFlag<0) return(-1);
  if (m->used==0) return(-2);

  /* if message not ready, deschedule and return 0 */
  if (!m->ready)
  {
    TimeWait(TimeNow()+DELAY);
    return(0);
  }

  /* ok ready */
  m->used = 0;
  return(1);
}


VChannelPtr ConnASync (int p, int id)
{
  int i,error;
  Slot_t *mySlot;

  /* find free slot */
  Wait(&mySema);
  if (usedCount<MAXSLOTS)
  {
    for (i=0; i<MAXSLOTS; i++)
      if (Slots[slotPtr].used)
        slotPtr = (slotPtr+1)%MAXSLOTS;
      else
        break;
    if (i>=MAXSLOTS)
    {
      Signal(&mySema);
      return(NULL);
    }
    mySlot = &(Slots[slotPtr]);
    mySlot->used = 1;
    usedCount++;
    slotPtr = (slotPtr+1)%MAXSLOTS;
    Signal(&mySema);
  }
  else
  {
    Signal(&mySema);
    printe("%d: no more slots in ConnectASync(), dest=%d, id=%d\n",me,p,id);
    return(NULL);
  }

  mySlot->SendFlag = 1;
  mySlot->RecvFlag = 1;

  /* init the two fifos */
  error = Fifo_init(&(mySlot->outgoing));
  if (error!=0)
  {
    printe("%d: error in ConnectASync() initializing outgoing fifo, dest=%d, id=%d\n",me,p,id);
    mySlot->used = 0;
    return(NULL);
  }
  error = Fifo_init(&(mySlot->incoming));
  if (error!=0)
  {
    printe("%d: error in ConnectASync() initializing incoming fifo, dest=%d, id=%d\n",me,p,id);
    Fifo_exit(&(mySlot->outgoing));
    mySlot->used = 0;
    return(NULL);
  }

  /* start sender thread */
  mySlot->SendThread = StartThread((FunctionPtr)SenderThread,STACKSIZE,&error,
                                   2*sizeof(int)+sizeof(mySlot),p,id,mySlot);
  if (mySlot->SendThread==NULL)
  {
    printe("%d: could not start SenderThread in ConnectASync(), dest=%d, id=%d, error = %d\n",me,p,id,error);
    mySlot->SendFlag = -99;
  }

  /* start receiver thread */
  mySlot->RecvThread = StartThread((FunctionPtr)ReceiverThread,STACKSIZE,&error,
                                   3*sizeof(int)+sizeof(mySlot),p,id,mySlot);
  if (mySlot->RecvThread==NULL)
  {
    printe("%d: could not start ReceiverThread in ConnectASync(), dest=%d, id=%d, error = %d\n",me,p,id,error);
    mySlot->RecvFlag = -99;
  }

  /* return ok */
  return((VChannelPtr) mySlot);
}

int InfoAConn (VChannelPtr vc)
{
  int result,sf,rf;
  Slot_t *mySlot;

  /* cast to slot pointer */
  mySlot = (Slot_t *) vc;

  sf = mySlot->SendFlag;
  rf = mySlot->RecvFlag;

  if (VERBOSE>=2)
  {
    printe("%d: InfoAConn: slot=%x, sf=%d, rf=%d\n",me,vc,sf,rf);
  }

  /* if both flags are positive, ok */
  if ((sf==0)&&(rf==0)) return(1);

  /* if one flag is 1, then still wait */
  if ((sf==1)||(rf==1))
  {
    TimeWait(TimeNow()+DELAY);
    return(0);
  }

  /* error condition, this is absolutely fatal */
  return(-1);
}


int DiscASync (VChannelPtr vc)
{
  MessageCarrier *messSend,*messRecv;
  Slot_t *mySlot;
  int i;

  /* cast to slot pointer */
  mySlot = (Slot_t *) vc;

  if (VERBOSE>=2)
  {
    printe("%d: DiscASync: slot=%x\n",me,mySlot);
  }

  /* allocate two message carriers */
  for (i=0; i<ASYNCMAX; i++)
    if (mc[mcPtr].used)
      mcPtr = (mcPtr+1)%ASYNCMAX;
    else
    {
      messSend = &(mc[mcPtr]);
      break;
    }
  for (i=0; i<ASYNCMAX; i++)
    if (mc[mcPtr].used)
      mcPtr = (mcPtr+1)%ASYNCMAX;
    else
    {
      messRecv = &(mc[mcPtr]);
      break;
    }
  if (i>=ASYNCMAX) return(1);
  messSend->used = 1;
  messRecv->used = 1;
  mcPtr = (mcPtr+1)%ASYNCMAX;

  /* clear flags */
  mySlot->SendFlag = 0;
  mySlot->RecvFlag = 0;

  /* send break messages */
  messSend->cmd = 1;
  messRecv->cmd = 1;
  Fifo_in(&(mySlot->outgoing),messSend);
  Fifo_in(&(mySlot->incoming),messRecv);

  /* return ok */
  return(0);
}

int InfoADisc (VChannelPtr vc)
{
  int sf,rf,sr,rr;
  Slot_t *mySlot;

  /* cast to slot pointer */
  mySlot = (Slot_t *) vc;

  sf = mySlot->SendFlag;
  rf = mySlot->RecvFlag;

  /* if one flag is zero, then still wait */
  if ((sf!=1)||(rf!=1))
  {
    TimeWait(TimeNow()+DELAY);
    return(0);
  }

  /* wait for both threads */
  WaitThread(mySlot->SendThread,&sr);
  WaitThread(mySlot->RecvThread,&rr);

  /* clear fifos, free slot */
  Fifo_exit(&(mySlot->outgoing));
  Fifo_exit(&(mySlot->incoming));
  Wait(&mySema);
  mySlot->used = 0;
  usedCount--;
  Signal(&mySema);

  if ((sr==0)&&(rr==0)) return(1);

  return(-1);
}

/****************************************************************************/
/*                                                                          */
/* Random communication                                                         */
/*                                                                          */
/****************************************************************************/

int SendMail (int destId, int reqId, void *data, int size)
{
  if (size>RAND_MSG_SIZE) return(1);
  return(-PutMessage(destId,reqId,MSGTYPE,1,-1,data,size));
}

int GetMailProcess (void)
{
  RR_Message_t Message;
  RandomMessage *rmslot;
  int i;

  ChangePriority(LOW_PRIORITY);
  while (1)
  {
    i = GetMessage(-1,-1,MSGTYPE,-1,&Message);
    if (rmExit) return(0);
    Wait(&rmSema);
    if (i<0)
    {
      rmError++;
      Signal(&rmSema);
      continue;
    }
    if (Message.Header.Size>RAND_MSG_SIZE)
    {
      rmError++;
      Signal(&rmSema);
      continue;
    }
    if (rmCount<MAILBOXSIZE)
    {
      for (i=0; i<MAILBOXSIZE; i++)
        if (rm[rmRecvPtr].used)
          rmRecvPtr = (rmRecvPtr+1)%MAILBOXSIZE;
        else
          break;
      if (i>=MAILBOXSIZE)
      {
        rmError++;
        Signal(&rmSema);
        continue;
      }
      rmslot = &(rm[rmRecvPtr]);
      rmslot->used = 1;
      rmCount++;
      rmRecvPtr = (rmRecvPtr+1)%MAILBOXSIZE;
      Signal(&rmSema);
      rmslot->sourceId = Message.Header.SourceProcId;
      rmslot->reqId = Message.Header.ReqId;
      rmslot->size = Message.Header.Size;
      memcpy(rmslot->data,Message.Body,rmslot->size);
      continue;
    }
    rmError++;
    Signal(&rmSema);
  }
  return(0);
}

int GetMail (int *sourceId, int *reqId, void *data, int *size)
{
  int i;
  RandomMessage *rmslot;

  Wait(&rmSema);

  /* if errors occured */
  if (rmError>0)
  {
    i = rmError;
    rmError = 0;
    Signal(&rmSema);
    return(-i);                                 /* return neg. value */
  }

  /* if there are messages in the box */
  if (rmCount>0)
  {
    for (i=0; i<MAILBOXSIZE; i++)
      if (rm[rmReadPtr].used)
        break;
      else
        rmReadPtr = (rmReadPtr+1)%MAILBOXSIZE;
    if (i>=MAILBOXSIZE)
    {
      Signal(&rmSema);
      return(-1);
    }
    rmslot = &(rm[rmReadPtr]);
    rmslot->used = 0;
    i = rmCount--;
    rmReadPtr = (rmReadPtr+1)%MAILBOXSIZE;
    *sourceId = rmslot->sourceId;
    *reqId = rmslot->reqId;
    *size = rmslot->size;
    memcpy(data,rmslot->data,rmslot->size);
    Signal(&rmSema);
    return(i);
  }

  /* no messages in the box */
  Signal(&rmSema);
  return(0);
}


/****************************************************************************/
/*                                                                          */
/* Miscellaneous                                                                        */
/*                                                                          */
/****************************************************************************/

int UsedSpace (void)
{
  return((int)(100.0*((float)usedCount)/((float)MAXSLOTS)));
}

void PrintHostMessage (char *s)
{
  printf("%s",s);
}

double CurrentTime (void)
{
  return( ((double)TimeNow())/((double)CLOCK_TICK) );
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
