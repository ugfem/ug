// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      pvm-ppif.c                                                    */
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
/*            PVM 3  Implementation                                         */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*            internet: bastian@iwr1.iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   13 Sep 1993, begin                                            */
/*            951106 kb  changed parameters of InitPPIF()                   */
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
#include <time.h>

#include "pvm3.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define RECVMAX     128         /* max number of recv messages at one time      */
#define MAXPROC     128     /* max number of processes allowed              */
#define VERBOSE         0               /* verbose level								*/
#define DELAY           5000    /* delay in async receive						*/
#define PVMMODE     PvmDataRaw    /* set to PvmDataInPlace when available   */

#define MAXT        15      /* maximum number of downtree nodes max log2(P) */
#define ARRAYID     10000   /* channel id's for 3D array structure          */
#define TREEID      10001   /* channel id's for tree structure              */
#define RANDOMID    20000   /* msgtag offset for random messages            */

#define ABS(i)      (((i)<0) ? (-(i)) : (i))
#define MAKEVC(p,i) (((p)<<16)|(((i)+1)&0xFFFF))
#define GETPROC(vc) ((((vc)&0xFFFF0000)>>16)&0xFFFF)
#define GETID(vc)   (((vc)&0x0FFFF)-1)
#define XPOS(aid)   (aid&0xFF)                   /* xpos from compact form  */
#define YPOS(aid)   ((aid&0xFF00)>>8)            /* ypos from compact form  */
#define ZPOS(aid)   ((aid&0xFF0000)>>16)         /* zpos from compact form  */
#define MIN(x,y) (((x)<(y)) ? (x) : (y))

/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

typedef unsigned int VChannelPtr;
typedef struct message *msgid;

enum directions {north,east,south,west,up,down};

struct message {
  struct message *next,*prev;
  VChannelPtr vc;
  int size;
  void *data;
}  ;

typedef struct message MESSAGE ;


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

static MESSAGE pool[RECVMAX];     /* stores message parameters				*/
static MESSAGE *firstM=NULL;
static MESSAGE *lastM=NULL;       /* list of messages not yet received      */
static MESSAGE *freeM=NULL;       /* free list                              */

static int mytid;                 /* my task identifier                     */
static int tids[MAXPROC];         /* task identifiers                       */

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

int SendSync (VChannelPtr vc, void *data, int size);
int RecvSync (VChannelPtr vc, void *data, int size);

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
    case 1 :
      downtree[degree] = MAKEVC(aid_to_pid(MyX+1,MyY,MyZ),TREEID); break;
    case 2 :
      downtree[degree] = MAKEVC(aid_to_pid(MyX,MyY+1,MyZ),TREEID); break;
    case 3 :
      downtree[degree] = MAKEVC(aid_to_pid(MyX,MyY,MyZ+1),TREEID); break;
    default : return(1);
    }
    degree++;
  }

  /* uptree connections */
  if (((len==2)&&(mypos==start+1))||((len==3)&&(mypos>=start+1)))
  {
    switch (dim)
    {
    case 1 : uptree = MAKEVC(aid_to_pid(MyX-1,MyY,MyZ),TREEID); break;
    case 2 : uptree = MAKEVC(aid_to_pid(MyX,MyY-1,MyZ),TREEID); break;
    case 3 : uptree = MAKEVC(aid_to_pid(MyX,MyY,MyZ-1),TREEID); break;
    default : return(1);
    }
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
        case 1 : downtree[degree] = MAKEVC(aid_to_pid(start+l1,MyY,MyZ),TREEID); break;
        case 2 : downtree[degree] = MAKEVC(aid_to_pid(MyX,start+l1,MyZ),TREEID); break;
        case 3 : downtree[degree] = MAKEVC(aid_to_pid(MyX,MyY,start+l1),TREEID); break;
        default : return(1);
        }
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
        case 1 : uptree = MAKEVC(aid_to_pid(start,MyY,MyZ),TREEID); break;
        case 2 : uptree = MAKEVC(aid_to_pid(MyX,start,MyZ),TREEID); break;
        case 3 : uptree = MAKEVC(aid_to_pid(MyX,MyY,start),TREEID); break;
        default : return(1);
        }
      }
      return(BuildTree(start+l1,l2,dim));
    }
  }

  /* return ok */
  return(0);
}

static void fit_configuration (int p, int *xx, int *yy)
{
  int i,x,y;

  x = y = 1;

  for (i=0; i<p; i++)
  {
    if ((x<=y)&&((x+1)*y<=p))
    {
      x++;
      continue;
    }
    if ((y<x)&&(x*(y+1)<=p))
    {
      y++;
      continue;
    }
  }
  *xx = x; *yy = y;
  return;
}

int InitPPIF (int *argcp, char ***argvp)
{
  int i,j,error,cnt,configok;
  int nhost,narch;
  struct pvmhostinfo *hostp;

  /* get my id */
  mytid = pvm_mytid();
  if (mytid<0) return(1);
  /* pvm_serror(1); */

  /* find out if parent or child */
  tids[0] = pvm_parent();
  if (tids[0]<0)
  {
    tids[0] = mytid;
    me = 0;

    /* find configuration */
    configok = 0;
    DimZ = 1;
    for (i=1; i<*argcp; i++)
      if (strncmp((*argvp)[i],"-sz",2)==0)
      {
        int j;

        if (i+2>=*argcp)
        {
          printf("not enough arguments for -sz option\n");
          break;
        }
        sscanf((*argvp)[i+1],"%d",&DimX);
        sscanf((*argvp)[i+2],"%d",&DimY);
        if ((DimX<1)||(DimY<1)||(DimX*DimY>MAXPROC)) break;
        configok = 1;

        /* remove used arguments from arglist, 960409 KB */
        *argcp = *argcp - 3;
        for(j=i; j<*argcp; j++)
          (*argvp)[j] = (*argvp)[j+3];
        break;
      }
    if (!configok)
    {
      error = pvm_config(&nhost,&narch,&hostp);
      if (error<0)
      {
        printf("Error in reading configuration\n");
        return(1);
      }
      if (narch>1)
      {
        printf("Sorry, please use only one type of machine\n");
        return(1);
      }
      fit_configuration(nhost,&DimX,&DimY);
    }

    /* spawn processes */
    procs = DimX*DimY;
    printf("Spawning %s on %d by %d processors\n",(*argvp)[0],DimX,DimY);

    /* now start ugp on other processors */
    if (procs>1)
    {
      error=pvm_spawn((*argvp)[0],(char **) 0,PvmTaskDefault,"",procs-1,&tids[1]);
      if (error<procs-1)
      {
        printf("Error in spawning tasks error=%d\n",error);
        for (i=0; i<procs; i++)
          printf("tids[%d]=%d\n",i,tids[i]);
        return(1);
      }

      /* send DimX */
      pvm_initsend( PVMMODE );
      pvm_pkint(&DimX,1,1);
      pvm_mcast(&tids[1],procs-1,0);

      /* send DimY */
      pvm_initsend( PVMMODE );
      pvm_pkint(&DimY,1,1);
      pvm_mcast(&tids[1],procs-1,1);

      /* send tid array */
      pvm_initsend( PVMMODE );
      pvm_pkint(tids,procs,1);
      pvm_mcast(&tids[1],procs-1,2);
    }
  }
  else
  {
    /* receive DimX */
    pvm_recv(tids[0],0);
    pvm_upkint(&DimX,1,1);

    /* receive DimY */
    pvm_recv(tids[0],1);
    pvm_upkint(&DimY,1,1);

    /* set Dim and procs */
    DimZ = 1;
    procs = DimX*DimY;

    /* receive tid array */
    pvm_recv(tids[0],2);
    pvm_upkint(tids,procs,1);

    /* find my position */
    for (i=1; i<procs; i++)
      if (mytid==tids[i]) {me = i; break;}
  }


  /* compute additional config parameters */
  master = 0;
  arrayid = pid_to_aid(me);
  MyX = XPOS(arrayid);
  MyY = YPOS(arrayid);
  MyZ = ZPOS(arrayid);

  /* initialize array topology */
  for (i=0; i<6; i++) nn[i] = (VChannelPtr) NULL;
  if (DimX>1)
  {
    if (MyX>0) nn[west] = MAKEVC(me-1,ARRAYID);
    if (MyX<DimX-1) nn[east] = MAKEVC(me+1,ARRAYID);
  }
  if (DimY>1)
  {
    if (MyY>0) nn[south] = MAKEVC(me-DimX,ARRAYID);
    if (MyY<DimY-1) nn[north] = MAKEVC(me+DimX,ARRAYID);
  }
  if (DimZ>1)
  {
    if (MyZ>0) nn[down] = MAKEVC(me-DimX*DimY,ARRAYID);
    if (MyZ<DimZ-1) nn[up] = MAKEVC(me+DimX*DimY,ARRAYID);
  }

  /* tree configuration */
  degree = 0;
  uptree = 0;
  for (i=0; i<MAXT; i++) downtree[i] = 0;
  /* allocate uptree channel */
  if ((MyX==0)&&(MyY==0))
  {
    error = BuildTree(0,DimZ,3);
    if (error!=0)
    {
      printf("%d: error in BuildTree z-direction\n",me);
      return(1);
    }
  }
  if (MyX==0)
  {
    error = BuildTree(0,DimY,2);
    if (error!=0)
    {
      printf("%d: error in BuildTree y-direction\n",me);
      return(1);
    }
  }
  error = BuildTree(0,DimX,1);
  if (error!=0)
  {
    printf("%d: error in BuildTree x-direction\n",me);
    return(1);
  }

  /* count subtree nodes */
  for (i=0; i<MAXT; i++) slvcnt[i] = 0;
  for (i=0; i<degree; i++)
  {
    if (RecvSync(downtree[i],&j,sizeof(int))<0) return(1);
    slvcnt[i] = j;
  }
  cnt = 1;
  for (i=0; i<degree; i++) cnt += slvcnt[i];
  if (me!=master)
    if (SendSync(uptree,&cnt,sizeof(int))<0) return(1);

  /* init message pool */
  freeM = &(pool[0]);
  for (i=0; i<RECVMAX-1; i++) pool[i].next = &(pool[i+1]);
  pool[RECVMAX-1].next = NULL;

  return(0);
}

void ExitPPIF (void)
{
  pvm_exit();
  return;
}

/****************************************************************************/
/*                                                                          */
/* Tree oriented functions                                                  */
/*                                                                          */
/****************************************************************************/

int Broadcast (void *data, int size)
{
  int i,error,bufid;

  if (procs==1) return(0);

  if (me==master)
  {
    error = pvm_initsend(PVMMODE);
    if (error<0) return(-error);
    error = pvm_pkbyte(data,size,1);
    if (error<0) return(-error);
    if (VERBOSE)
    {
      printf("%4d: BCAST dst=xx id=%5d data=%8x size=%7d\n",me,TREEID,data,size);
    }
    error = pvm_mcast(&tids[1],procs-1,TREEID);
    if (error<0) return(-error);
  }
  else
  {
    if (VERBOSE)
    {
      printf("%4d: BCRCV src= 0 id=%5d data=%8x size=%7d\n",me,TREEID,data,size);
    }
    bufid = pvm_recv(tids[0],TREEID);
    if (bufid<0) return(1);
    pvm_upkbyte(data,size,1);
  }

  return(0);
}

int Concentrate (void *data, int size)
{
  if (me!=master)
    if (SendSync(uptree,data,size)<0) return(1);
  return(0);
}

int GetConcentrate (int slave, void *data, int size)
{
  if (slave<degree)
    if (RecvSync(downtree[slave],data,size)<0) return(1);
  return(0);
}


int Spread (int slave, void *data, int size)
{
  if (slave<degree)
  {
    if (SendSync(downtree[slave],data,size) != size) return(1);
  }
  return(0);
}

int GetSpread (void *data, int size)
{
  if (me!=master)
  {
    if (RecvSync(uptree,data,size) != size) return(1);
  }
  return(0);
}



int Synchronize (void)
{
  int i,n;

  /* send a dummy from leaves to root */
  for (i=degree-1; i>=0; i--) RecvSync(downtree[i],&n,sizeof(int));
  if (me!=master) SendSync(uptree,&n,sizeof(int));

  /* now all are there but only the root knows this */
  Broadcast(&n,sizeof(int));

  return(0);
}



/****************************************************************************/
/*                                                                          */
/* Synchronous communication                                                */
/*                                                                          */
/****************************************************************************/

VChannelPtr ConnSync (int p, int id)
{
  return(MAKEVC(p,id));
}

int DiscSync (VChannelPtr vc)
{
  return(0);
}

int SendSync (VChannelPtr vc, void *data, int size)
{
  int error;

  error = pvm_initsend(PVMMODE);
  if (error<0) return(-1);
  error = pvm_pkbyte(data,size,1);
  if (error<0) return(-1);
  if (VERBOSE)
  {
    printf("%4d: SSEND dst=%2d id=%5d data=%8x size=%7d\n",me,GETPROC(vc),GETID(vc),data,size);
  }
  error = pvm_send(tids[GETPROC(vc)],GETID(vc));
  if (error<0) return(-1);
  return(size);
}

int RecvSync (VChannelPtr vc, void *data, int size)
{
  int bufid,error;

  if (VERBOSE)
  {
    printf("%4d: SRECV src=%2d id=%5d data=%8x size=%7d\n",me,GETPROC(vc),GETID(vc),data,size);
  }
  bufid = pvm_recv(tids[GETPROC(vc)],GETID(vc));
  if (bufid<0) return(-1);
  error = pvm_upkbyte(data,size,1);
  if (error<0) return(-1);
  return(size);
}


/****************************************************************************/
/*                                                                          */
/* Asynchronous communication                                               */
/*                                                                          */
/****************************************************************************/

msgid SendASync (VChannelPtr vc, void *data, int size, int *error)
{
  int info;

  info = pvm_initsend(PVMMODE);
  if (info<0)
  {
    *error=1;
    printf("error in pvm_pkbyte vc=%x\n",vc);
    return(NULL);
  }
  info = pvm_pkbyte(data,size,1);
  if (info<0)
  {
    *error=1;
    printf("error in pvm_pkbyte vc=%x\n",vc);
    return(NULL);
  }
  if (VERBOSE)
  {
    printf("%4d: ASEND dst=%2d id=%5d data=%8x size=%7d\n",me,GETPROC(vc),GETID(vc),data,size);
  }
  info = pvm_send(tids[GETPROC(vc)],GETID(vc));
  if (info<0)
  {
    *error=1;
    printf("error in pvm_send vc=%x\n",vc);
    return(NULL);
  }
  *error = 0;

  return((msgid)42);
}

int InfoASend (VChannelPtr vc, msgid m)
{
  if (VERBOSE)
  {
    printf("%4d: ISEND dst=%2d id=%5d msid=%8x\n",me,GETPROC(vc),GETID(vc),m);
  }
  return(1);
}


msgid RecvASync (VChannelPtr vc, void *data, int size, int *error)
{
  MESSAGE *mess;

  /* allocate a free message carrier */
  mess = freeM;
  if (mess==NULL) {*error=1; return(NULL);}
  freeM = freeM->next;

  /* fill out MESSAGE structure */
  mess->vc = vc;
  mess->data = data;
  mess->size = size;

  /* append mess to list of messages to receive */
  if (lastM==NULL)
  {
    firstM = lastM = mess;
    mess->prev = mess->next = NULL;
  }
  else
  {
    lastM->next = mess;
    mess->prev = lastM;
    mess->next = NULL;
    lastM = mess;
  }

  *error = 0;
  if (VERBOSE)
  {
    printf("%4d: ARECV src=%2d id=%5d data=%8x size=%7d msid=%8x\n",me,GETPROC(vc),GETID(vc),data,size,mess);
  }

  return(mess);
}


int InfoARecv (VChannelPtr vc, msgid m)
{
  int bufid;
  MESSAGE *mess,*first;
  int info,bytes,msgtag,tid,found;
  int Case,i;

  /* is m the oldest message we are waiting for ? */
  found = 0; first=NULL;
  for (mess=firstM; mess!=NULL; mess=mess->next)
  {
    if ((mess->vc=vc)&&(first==NULL)) first = mess;
    if (mess==m) found++;
  }

  /* if no message for that channel exists,  fatal error */
  if (first==NULL)
  {
    printf("no message to receive from that channel vc=%x\n",vc);
    return(-1);
  }
  mess = first;

  if (found!=1)
  {
    printf("invalid msgid vc=%x msgid=%x\n",vc,m);
    return(-1);
  }

  /* check if m is the first message in the list */
  if (m!=first) {Case=1; goto contexit;}

  /* lets see if message has arrived */
  bufid = pvm_probe(tids[GETPROC(vc)],GETID(vc));
  if (bufid<0)
  {
    printf("pvm_probe failed in InfoARecv vc=%x msgid=%x\n",vc,m);
    return(-1);
  }
  if (bufid==0) {
    Case=2;
    goto contexit;
  }

  /* lets check size of messages */
  bufid = pvm_recv(tids[GETPROC(vc)],GETID(vc));
  if (bufid<0)
  {
    printf("pvm_recv failed in InfoARecv vc=%x msgid=%x\n",vc,m);
    return(-1);
  }
  info = pvm_bufinfo(bufid,&bytes,&msgtag,&tid);
  if (info<0)
  {
    printf("bufinfo failed info=%d\n",info);
    return(-1);
  }
  if (bytes!=mess->size)
  {
    printf("size inconsistent vc=%x bytes=%d size=%d\n",vc,bytes,mess->size);
    return(-1);
  }

  /* copy to ug's buffer */
  pvm_upkbyte(mess->data,mess->size,1);

  /* recycle MESSAGE structure */
  if (mess->next!=NULL) mess->next->prev = mess->prev;
  if (mess->prev!=NULL) mess->prev->next = mess->next;
  if (firstM==mess) firstM = mess->next;
  if (lastM==mess) lastM = mess->prev;
  mess->next = freeM;
  freeM = mess;
  Case=0;

  if (VERBOSE)
  {
    printf("%4d: IRECV src=%2d id=%5d msid=%8x stat=%1d\n",me,GETPROC(vc),GETID(vc),m,Case);
  }
  return(1);

contexit: /* continue polling */
  if (VERBOSE)
  {
    printf("%4d: IRECV src=%2d id=%5d msid=%8x stat=%1d\n",me,GETPROC(vc),GETID(vc),m,Case);
  }
  /* busy wait to prevent to many calls to the pvmd */
  i = 0;
  while (i<DELAY) i++;

  return(0);
}


VChannelPtr ConnASync (int p, int id)
{
  return(MAKEVC(p,id));
}

int InfoAConn (VChannelPtr vc)
{
  return(1);
}


int DiscASync (VChannelPtr vc)
{
  return(0);
}

int InfoADisc (VChannelPtr vc)
{
  return(1);
}

/****************************************************************************/
/*                                                                          */
/* Random communication                                                         */
/*                                                                          */
/****************************************************************************/

int SendMail (int destId, int reqId, void *data, int size)
{
  int error;
  error = pvm_initsend(PVMMODE);
  if (error<0) return(-1);
  error = pvm_pkbyte(data,size,1);
  if (error<0) return(-1);
  if (VERBOSE)
  {
    printf("%4d: SNDM  dst=%2d id=%5d data=%8x size=%7d\n",me,destId,reqId+RANDOMID,data,size);
  }
  error = pvm_send(tids[destId],reqId+RANDOMID);
  if (error<0)
  {
    printf("pvm_send failed in SendMail\n");
    return(-1);
  }

  return(0);
}

int GetMail (int *sourceId, int *reqId, void *data, int *size)
{
  int error,bufid,bytes,msgtag,tid,i;

  bufid = pvm_probe(-1,-1);
  if (bufid==0) return(0);
  error = pvm_bufinfo(bufid,&bytes,&msgtag,&tid);
  if (error<0)
  {
    printf("pvm_bufinfo failed in GetMail\n");
    return(-1);
  }
  if (msgtag<RANDOMID) return(0);

  /* receive message */
  for (i=0; i<procs; i++) if (tids[i]==tid) {*sourceId=i; break;}
  if (VERBOSE)
  {
    printf("%4d: GETM  src=%2d id=%5d data=%8x size=%d\n",me,*sourceId,msgtag,data,bytes);
  }
  error = pvm_recv(tid,msgtag);
  if (error<0)
  {
    printf("pvm_recv failed in GetMail\n");
    return(-1);
  }
  error = pvm_upkbyte(data,bytes,1);
  if (error<0)
  {
    printf("pvm_upkbyte failed in GetMail\n");
    return(-1);
  }

  /* fill parameters */
  *size = bytes;
  *reqId = msgtag-RANDOMID;

  return(1);
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

float CurrentTime (void)
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

  return(ABS(pX-qX)+ABS(pY-qY)+ABS(pZ-qZ));
}
