// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      T3D/ppif.c                                                    */
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
/*            CRAY T3D Implementation                                       */
/*                                                                          */
/* Author:    Stefan Lang                                                   */
/*            Institut fuer Computeranwendungen 3 ( ICA3 )                  */
/*            Pfaffenwaldring 27                                            */
/*            Universitaet Stuttgart                                        */
/*            70569 Stuttgart                                               */
/*            e-mail: stefan@ica3.uni-stuttgart.de                          */
/*                                                                          */
/* History:   8 Dec 1994, begin                                             */
/*                                                                          */
/* Remarks:                                                                 */
/*            -number of PEs should be a multiple of 2                      */
/*            -T3D architecture requires to specify a even number of PEs    */
/*             in x dimension                                               */
/*            -current configuration is 4x4x2                               */
/*            - currently allowed convex shapes are:                        */
/*                                                                          */
/*                          PEs    x   y   z                                */
/*                                                                          */
/*                          1(?)  1   1   1                                 */
/*                          2     2   1   1                                 */
/*                          4     2   2   1                                 */
/*                                2   1   2                                 */
/*                          6     2   3   1                                 */
/*                          8     2   2   2                                 */
/*                                4   2   1                                 */
/*                                4   1   2                                 */
/*                                2   4   1                                 */
/*                          12    4   3   1                                 */
/*                          16    4   2   2                                 */
/*                                4   4   1                                 */
/*                                2   4   2                                 */
/*                          24    4   3   2                                 */
/*                          32    4   4   2                                 */
/*                                                                          */
/*             (for more information see Cray MPP User's Guide SG-2508 1.1) */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <mpp/shmem.h>
#include <mpp/sync_proto.h>

#include "../ppif_general.h"

#include "compiler.h"


/* TODO: delete this */
/* copied from compiler.h */
/*#define ALIGNMENT 8   */                  /* power of 2 and >= sizeof(int) !  */
/*#define ALIGNMASK 0xFFFFFFF8  */          /* compatible to alignment          */


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/* macros for conditional compilation */
#define PCKSIZE4BYTE          /* data is transmitted in 4 byte packagesi, if  */
                              /* 8 byte pckgs wanted comment this define out  */
#ifdef PCKSIZE4BYTE
#define PACKSIZE          4   /* value must be depend of PCKSIZE4BYTE 4 or 8  */
#else
#define PACKSIZE          8   /* value must be depend of PCKSIZE4BYTE 4 or 8  */
#endif

#define MAXPROC         256   /* max number of processors                     */
#define MAXT              8   /* maximum number of downtree nodes max log2(P) */
#define SYNCMAX         256   /* maximum number of sync virtual channels      */
#define ASYNCMAX        256   /* maximum number of async virtual channels     */
#define ASYNCSENDMAX     64   /* max of send messages at one time             */
#define ASYNCRECEIVEMAX  64   /* max of receive messages at one time          */
#define MAILMAX         256   /* max of mails one PE cna handle at one time   */
#define RAND_MSG_SIZE    32   /* max size in bytes of one mail                */
#define ARRAYID     10000     /* channel id's for 3D array structure          */
#define TREEID      10001     /* channel id's for tree structure              */

/*#define TIME  */                /* timing level (0 = no timing)                 */
#define VERBOSE           0   /* verbose level for general output             */
#define SYNCVERBOSE       0   /* verbose level for sync com                   */
#define ASYNCVERBOSE      0   /* verbose level for async com                  */
#define NOPE             -1   /* value if channel entry is not used           */
#define NOID             -1   /* value if channel id is undefined             */
#define NOINDEX          -1   /* value if channel has no corresponding entry  */
#define NOIDENT          -1   /* value if ident is undefined                  */
#define NOSIZE           -1   /* value if size is undefined                   */
#define RECEIVED         -1   /* value receiver has received                  */
#define NOTRECEIVED       1   /* value receiver has not received              */
#define OCCUPIED          1   /* value if mail entry is occupied              */
#define NOTOCCUPIED       0   /* value if mail entry is not occupied = free   */
#define READABLE          1   /* value if mail entry is readable              */
#define NOTREADABLE       0   /* value if mail entry is not readable          */
#define NOFIRST          -1   /* value if first is still searched             */
#define NODATA            0   /* value if data array is empty                 */
#define NOFLAG           -1   /* value if flag entry id undefined             */

/* general macros */
#define XPOS(aid)   (aid&0xFF)                   /* xpos from compact form  */
#define YPOS(aid)   ((aid&0xFF00)>>8)            /* ypos from compact form  */
#define ZPOS(aid)   ((aid&0xFF0000)>>16)         /* zpos from compact form  */
#define ABS(i)      (((i)<0) ? (-(i)) : (i))


/* macros for sync comm */
#define SYNCVCHANNELSIZE                 (SYNCMAX*sizeof(SYNCVCHANNEL))
#define SET_SVC_PE(svc,p)                ((svc)->pe = (p))
#define SET_SVC_ID(svc,i)                ((svc)->id = (i))
#define SET_SVC_INDEX(svc,j)             ((svc)->index = (j))
#define SET_SVC_RECEIVED(svc,j)          ((svc)->received = (j))
#define SET_SVC_DATA(svc,j)              ((svc)->data = (j))
#define SET_SVC_SIZE(svc,j)              ((svc)->size = (j))
#define GET_PE(svc)                      ((svc)->pe)
#define GET_ID(svc)                      ((svc)->id)
#define GET_INDEX(svc)                   ((svc)->index)
#define GET_SIZE(svc)                    ((svc)->size)
#define GET_SVCA_ADR(sync_array,i)       (&(sync_array)[(i)])
#define GET_SVCA_PE(i)                   (sync_channels[(i)].pe)
#define GET_SVC_DATA(svc)                ((svc)->data)
#define GET_SVC_DATAADR(svc)             (&sync_channels[(svc)->index].data)
#define GET_SVC_SIZEADR(svc)             (&sync_channels[(svc)->index].size)
#define GET_SVC_SRECVADR(svc)            (&sync_channels[(svc)->index].received)
#define GET_SVC_RECVADR(svc)             (&(svc)->received)
#define GET_SVC_ID(sync_channels,i)      ((sync_channels)[(i)].id)
#define MAKE_IDENT(id,index)             (((id)<<16)|((index)&0xFFFF))
#define GET_IDENT_ID(ident)              ((((ident)&0xFFFF0000)>>16)&0xFFFF)
#define GET_IDENT_INDEX(ident)           ((ident)&0x0000FFFF)

/* macros for async comm */
#define ASYNCVCHANNELSIZE                (ASYNCMAX*sizeof(ASYNCVCHANNEL))
#define SET_AVC_PE(avc,p)                ((avc)->pe = (p))
#define SET_AVC_ID(avc,i)                ((avc)->id = (i))
#define SET_AVC_INDEX(avc,j)             ((avc)->index = (j))
#define SET_AVC_NEXTSEND(avc,j)          ((avc)->next_send = (j))
#define SET_AVC_NEXTRECEIVE(avc,j)       ((avc)->next_receive = (j))
#define SET_AVC_CARRIERS(avc,j)          ((avc)->carriers = &async_carriers[(j)])
#define SET_CARRIER_DATA(car,ptr)        ((car)->data = (ptr))
#define SET_CARRIER_SIZE(car,i)          ((car)->size = (i))
#define SET_CARRIER_RECEIVED(car,i)      ((car)->received = (i))
#define GET_AVC_ADR(i)                   (&async_channels[(i)])
#define GET_AVC_PE(avca,i)                       ((avca)[(i)].pe)
#define GET_AVC_PEADR(avca,i)                    (&(avca)[(i)].pe)
#define GET_AVC_ID(avca,i)               ((avca)[(i)].id)
#define GET_AVC_IDADR(avca,i)                    (&(avca)[(i)].id)
#define GET_AVC_INDEX(avc)               (avc->index)
#define GET_AVC_INDEXADR(avcarray,j)     (&(avcarray)[(j)].index)
#define GET_AVC_CNTSEND(avc)             ((avc)->cnt_send)
#define GET_AVC_CNTRECEIVE(avc)          ((avc)->cnt_receive)
#define GET_AVC_NEXTSEND(avc)            ((avc)->next_send)
#define GET_AVC_NEXTRECEIVE(avc)         ((avc)->next_receive)
#define GET_AVC_CARRIERS(avc)            ((avc)->carriers)
#define GET_SCARRIER_ADR(avc,index)      (&(avc)->carriers->out_carriers[(index)])
#define GET_RCARRIER_ADR(avc,index)      (&(avc)->carriers->in_carriers[(index)])
#define GET_CARRIER_DATA(car)            ((car)->data)
#define GET_CARRIER_SIZE(car)            ((car)->size)
#define GET_CARRIER_RECEIVED(car)        ((car)->received)

/* macros for random communication */
#define MAILSIZE                         (sizeof(MAIL))
#define NEXT_MAIL_INDEX(i)               (((i)+1)%MAILMAX)
#define SET_MAIL_OCCUPIED(mail_ptr,i)    ((mail_ptr)->occupied = (i))
#define SET_MAIL_SOURCE(mail_ptr,i)      ((mail_ptr)->source = (i))
#define SET_MAIL_ID(mail_ptr,i)          ((mail_ptr)->id = (i))
#define SET_MAIL_SIZE(mail_ptr,i)        ((mail_ptr)->size = (i))
#define SET_MAIL_READABLE(mail_ptr,i)    ((mail_ptr)->readable = (i))
#define GET_MAIL_READABLE(i)             (mailbox.mailarray[(i)].readable)
#define GET_MAIL_ADRESS(first)           (&mailbox.mailarray[(first)])
#define GET_MAIL_SOURCE(mail_ptr)        ((mail_ptr)->source)
#define GET_MAIL_ID(mail_ptr)            ((mail_ptr)->id)
#define GET_MAIL_SIZE(mail_ptr)          ((mail_ptr)->size)
#define GET_MAIL_DATA(mail_ptr)          ((mail_ptr)->data)
#define GET_MAIL_OCCUPIED(mail_ptr)      ((mail_ptr)->occupied)


/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/

struct SyncVChannel {
  int pe;
  int id;
  int index;
  int received;
  void *data;
  int size;
} ;

typedef struct SyncVChannel SYNCVCHANNEL;

struct carrier {
  int received;
  int size;
  void *data;
} ;

typedef struct carrier CARRIER;

typedef CARRIER *msgid;

struct ASyncVCCarrier {
  CARRIER out_carriers[ASYNCSENDMAX];
  CARRIER in_carriers[ASYNCRECEIVEMAX];
} ;

typedef struct ASyncVCCarrier ASYNCVCCARRIER;

struct ASyncVChannel {
  int pe;
  int id;
  int index;
  int next_send;
  int next_receive;
  ASYNCVCCARRIER *carriers;
} ;

typedef struct ASyncVChannel ASYNCVCHANNEL;

typedef int *VChannelPtr;

struct Mail {
  int occupied;
  int source;
  int id;
  int size;
  char data[RAND_MSG_SIZE];
  int readable;
} ;

typedef struct Mail MAIL;

struct Mailbox {
  int first_read;
  int first_write;
  MAIL mailarray[MAILMAX];
} ;

typedef struct Mailbox MAILBOX;

enum directions {north,east,south,west,up,down};


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
int arrayid;                /* compact format of position, 8 bits each      */
int MyX,MyY,MyZ;            /* 3D array coordinates                         */
int DimX,DimY,DimZ;         /* 3D array dimensions, may be 1 !              */
SYNCVCHANNEL nnsvca[6];     /* nearest neighbor channels in 3D array        */
SYNCVCHANNEL *nn[6];        /* pointers to nearest neighbor channels        */

/* Tree structure */
int degree;                   /* degree of downtree nodes                     */
VChannelPtr uptree;           /* channel uptree                               */
VChannelPtr downtree[MAXT];   /* channels downtree (may be empty)             */
int slvcnt[MAXT];             /* number of processors in subtree              */

#ifdef TIME
double timebegin,timeend;
#endif

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* definitions for channels of tree structure */
static SYNCVCHANNEL treesvca[MAXT+1];  /* channels downtree (may be empty)  */

/* definitions for synchronous communication */
static int sync_flag[MAXPROC];                    /* flags for synchronisation */
static SYNCVCHANNEL sync_channels[SYNCMAX+MAXT+1]; /* array of sync  VChannels  */

/* definitions for asynchronous communication */
static ASYNCVCHANNEL async_channels[ASYNCMAX];
static ASYNCVCHANNEL remote_channels[ASYNCMAX];
static ASYNCVCCARRIER async_carriers[ASYNCMAX];
static CARRIER carrier[1];

/* definitions for random communication */
static MAILBOX mailbox;
static MAIL newmail[1];

/* definition of clocks per second */
long clocks;

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

int SendSync (VChannelPtr vc, void *data, int size);
int RecvSync (VChannelPtr vc, void *data, int size);
double CurrentTime (void);


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
  SYNCVCHANNEL *adress;

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
    case 1 : adress = GET_SVCA_ADR(sync_channels,degree+1);
      downtree[degree] = (VChannelPtr) adress;
      SET_SVC_PE(adress,aid_to_pid(MyX+1,MyY,MyZ));
      SET_SVC_ID(adress,TREEID);
      SET_SVC_INDEX(adress,0);
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
      shmem_swap((long *) &GET_INDEX(GET_SVCA_ADR(sync_channels,0)),degree+1,GET_PE(adress));
      break;
    case 2 : adress = GET_SVCA_ADR(sync_channels,degree+1);
      downtree[degree] = (VChannelPtr) adress;
      SET_SVC_PE(adress,aid_to_pid(MyX,MyY+1,MyZ));
      SET_SVC_ID(adress,TREEID);
      SET_SVC_INDEX(adress,0);
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
      shmem_swap((long *) &GET_INDEX(GET_SVCA_ADR(sync_channels,0)),degree+1,GET_PE(adress));
      break;
    case 3 : adress = GET_SVCA_ADR(sync_channels,degree+1);
      downtree[degree] = (VChannelPtr) adress;
      SET_SVC_PE(adress,aid_to_pid(MyX,MyY,MyZ+1));
      SET_SVC_ID(adress,TREEID);
      SET_SVC_INDEX(adress,0);
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
      shmem_swap((long *) &GET_INDEX(GET_SVCA_ADR(sync_channels,0)),degree+1,GET_PE(adress));
      break;
    default : return(1);
    }
    degree++;
  }

  /* uptree connections */
  if (((len==2)&&(mypos==start+1))||((len==3)&&(mypos>=start+1)))
  {
    switch (dim)
    {
    case 1 : adress = GET_SVCA_ADR(sync_channels,0);
      uptree = (VChannelPtr) adress;
      SET_SVC_PE(adress,aid_to_pid(MyX-1,MyY,MyZ));
      SET_SVC_ID(adress,TREEID);
      /* index is set by father */
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
      break;
    case 2 : adress = GET_SVCA_ADR(sync_channels,0);
      uptree = (VChannelPtr) adress;
      SET_SVC_PE(adress,aid_to_pid(MyX,MyY-1,MyZ));
      SET_SVC_ID(adress,TREEID);
      /* index is set by father */
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
      break;
    case 3 : adress = GET_SVCA_ADR(sync_channels,0);
      uptree = (VChannelPtr) adress;
      SET_SVC_PE(adress,aid_to_pid(MyX,MyY,MyZ-1));
      SET_SVC_ID(adress,TREEID);
      /* index is set by father */
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
      break;
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
        case 1 : adress = GET_SVCA_ADR(sync_channels,degree+1);
          downtree[degree] = (VChannelPtr) adress;
          SET_SVC_PE(adress,aid_to_pid(start+l1,MyY,MyZ));
          SET_SVC_ID(adress,TREEID);
          SET_SVC_INDEX(adress,0);
          SET_SVC_RECEIVED(adress,NOTRECEIVED);
          SET_SVC_DATA(adress,NULL);
          SET_SVC_SIZE(adress,NOSIZE);
          shmem_swap((long *) &GET_INDEX(GET_SVCA_ADR(sync_channels,0)),degree+1,GET_PE(adress));
          break;
        case 2 : adress = GET_SVCA_ADR(sync_channels,degree+1);
          downtree[degree] = (VChannelPtr) adress;
          SET_SVC_PE(adress,aid_to_pid(MyX,start+l1,MyZ));
          SET_SVC_ID(adress,TREEID);
          SET_SVC_INDEX(adress,0);
          SET_SVC_RECEIVED(adress,NOTRECEIVED);
          SET_SVC_DATA(adress,NULL);
          SET_SVC_SIZE(adress,NOSIZE);
          shmem_swap((long *) &GET_INDEX(GET_SVCA_ADR(sync_channels,0)),degree+1,GET_PE(adress));
          break;
        case 3 : adress = GET_SVCA_ADR(sync_channels,degree+1);
          downtree[degree] = (VChannelPtr) adress;
          SET_SVC_PE(adress,aid_to_pid(MyX,start+l1,MyZ));
          SET_SVC_ID(adress,TREEID);
          SET_SVC_INDEX(adress,0);
          SET_SVC_RECEIVED(adress,NOTRECEIVED);
          SET_SVC_DATA(adress,NULL);
          SET_SVC_SIZE(adress,NOSIZE);
          shmem_swap((long *) &GET_INDEX(GET_SVCA_ADR(sync_channels,0)),degree+1,GET_PE(adress));
          break;
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
        case 1 : adress = GET_SVCA_ADR(sync_channels,0);
          uptree = (VChannelPtr) adress;
          SET_SVC_PE(adress,aid_to_pid(start,MyY,MyZ));
          SET_SVC_ID(adress,TREEID);
          /* index is set by father */
          SET_SVC_RECEIVED(adress,NOTRECEIVED);
          SET_SVC_DATA(adress,NULL);
          SET_SVC_SIZE(adress,NOSIZE);
          break;
        case 2 : adress = GET_SVCA_ADR(sync_channels,0);
          uptree = (VChannelPtr) adress;
          SET_SVC_PE(adress,aid_to_pid(MyX,start,MyZ));
          SET_SVC_ID(adress,TREEID);
          /* index is set by father */
          SET_SVC_RECEIVED(adress,NOTRECEIVED);
          SET_SVC_DATA(adress,NULL);
          SET_SVC_SIZE(adress,NOSIZE);
          break;
        case 3 : adress = GET_SVCA_ADR(sync_channels,0);
          uptree = (VChannelPtr) adress;
          SET_SVC_PE(adress,aid_to_pid(MyX,MyY,start));
          SET_SVC_ID(adress,TREEID);
          /* index is set by father */
          SET_SVC_RECEIVED(adress,NOTRECEIVED);
          SET_SVC_DATA(adress,NULL);
          SET_SVC_SIZE(adress,NOSIZE);
          break;
        default : return(1);
        }
      }
      return(BuildTree(start+l1,l2,dim));
    }
  }

  /* return ok */
  return(0);
}

static void fit_configuration (int p, int *x, int *y, int *z)
{
  /* determine dimensions with intrinsic function */
  PSHAPE(x,y,z);

  /* minimum partition size 2 PEs, fit for single PE */
  if (p==1) *x=*y=*z=1;
}

static void OLD_fit_configuration (int p, int *xx, int *yy, int *zz)
{
  int i,x,y,z;

  x = y = z = 1;

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
  *xx = x; *yy = y; *zz = z;
  return;
}

int InitPPIF (int *argcp, char ***argvp)
{
  int i,j,error,cnt,configok;
  SYNCVCHANNEL   *adress;
  ASYNCVCHANNEL   *asyncadress;
  ASYNCVCCARRIER *cararray;
  CARRIER        *carrier;
  MAIL           *mailadr;

  void *dummy1, *dummy2;

  /* clocks per second */
  clocks = sysconf(_SC_CLK_TCK);

  /* this is the fix proposed by cray */
  for (i=1; i<256; i++)
  {
    dummy1 = malloc(i*1024*1024);
    if (dummy1==NULL) {
      dummy1 = malloc((i-1)*1024*1024);
      dummy2 = malloc(32);
      free(dummy1);
      break;
    }
    free(dummy1);
  }

  if (dummy1!=NULL) free(dummy1);

  /* enable automatic cache invalidation */
  shmem_set_cache_inv();

  /* initialize basic config parameter */
  me = _my_pe();
  procs = _num_pes();
  master = 0;

  /* print master PE and number of procs */
  if (VERBOSE && me == master)
  {
    printf("master is PE %d,", me);
    printf("procs = %d\n", procs);
  }

  /* assemble command line options */
  for (i=1; i<*argcp; i++)
    if (strncmp((*argvp)[i],"-n",1)==0)
    {
      if (i+1>=*argcp)
      {
        printf("not enough arguments for -n option\n");
        break;
      }
      sscanf((*argvp)[i+1],"%d",&procs);
      if (procs>MAXPROC) break;
      configok = 1;
      break;
    }

  /* fit dimension parameters to virtual PE numbering */
  fit_configuration(procs,&DimX,&DimY,&DimZ);

  /* for now only two dimensions */
  DimY = DimY * DimZ;
  DimZ = 1;

  if (me == master)
    printf("Configuration: %d PEs in an %dx%dx%d array\n", procs, DimX, DimY, DimZ);

  /* compute additional config parameters */
  arrayid = pid_to_aid(me);
  MyX = XPOS(arrayid);
  MyY = YPOS(arrayid);
  MyZ = ZPOS(arrayid);

  /* initialize flag array */
  for (i=0; i<MAXPROC; i++) {
    sync_flag[i] = NOFLAG;
  }

  /* initialize synchronous vchannel array */
  for (i=0; i<SYNCMAX+MAXT+1; i++) {
    adress = GET_SVCA_ADR(sync_channels,i);
    SET_SVC_PE(adress,NOPE);
    SET_SVC_ID(adress,NOID);
    SET_SVC_INDEX(adress,NOINDEX);
    SET_SVC_RECEIVED(adress,NOTRECEIVED);
    SET_SVC_DATA(adress,NULL);
    SET_SVC_SIZE(adress,NOSIZE);
  }

  /* initialize array topology */
  for (i=0; i<6; i++) nn[i] = (long) NULL;
  if (DimX>1)
  {
    if (MyX>0)
    {
      adress = GET_SVCA_ADR(nnsvca,west);
      nn[west] = adress;
      SET_SVC_PE(adress,aid_to_pid(MyX-1,MyY,MyZ));
      SET_SVC_ID(adress,ARRAYID);
      SET_SVC_INDEX(adress,east);
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
    }
    else
    {
      adress = NULL;
      nn[west] = adress;
    }

    if (MyX<DimX-1)
    {
      adress = GET_SVCA_ADR(nnsvca,east);
      nn[east] = adress;
      SET_SVC_PE(adress,aid_to_pid(MyX+1,MyY,MyZ));
      SET_SVC_ID(adress,ARRAYID);
      SET_SVC_INDEX(adress,west);
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
    }
    else
    {
      adress = NULL;
      nn[west] = adress;
    }
  }
  if (DimY>1)
  {
    if (MyY>0)
    {
      adress = GET_SVCA_ADR(nnsvca,south);
      nn[south] = adress;
      SET_SVC_PE(adress,aid_to_pid(MyX,MyY-1,MyZ));
      SET_SVC_ID(adress,ARRAYID);
      SET_SVC_INDEX(adress,north);
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
    }
    else
    {
      adress = NULL;
      nn[west] = adress;
    }

    if (MyY<DimY-1)
    {
      adress = GET_SVCA_ADR(nnsvca,north);
      nn[north] = adress;
      SET_SVC_PE(adress,aid_to_pid(MyX,MyY+1,MyZ));
      SET_SVC_ID(adress,ARRAYID);
      SET_SVC_INDEX(adress,south);
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
    }
    else
    {
      adress = NULL;
      nn[west] = adress;
    }
  }
  if (DimZ>1)
  {
    if (MyZ>0)
    {
      adress = GET_SVCA_ADR(nnsvca,down);
      nn[down] = adress;
      SET_SVC_PE(adress,aid_to_pid(MyX,MyY,MyZ-1));
      SET_SVC_ID(adress,ARRAYID);
      SET_SVC_INDEX(adress,up);
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
    }
    else
    {
      adress = NULL;
      nn[west] = adress;
    }

    if (MyZ<DimZ-1)
    {
      adress = GET_SVCA_ADR(nnsvca,up);
      nn[up] = adress;
      SET_SVC_PE(adress,aid_to_pid(MyX,MyY,MyZ+1));
      SET_SVC_ID(adress,ARRAYID);
      SET_SVC_INDEX(adress,down);
      SET_SVC_RECEIVED(adress,NOTRECEIVED);
      SET_SVC_DATA(adress,NULL);
      SET_SVC_SIZE(adress,NOSIZE);
    }
    else
    {
      adress = NULL;
      nn[west] = adress;
    }
  }

  /* due to race conditions beginning at 256 procs up */
  /* Rudi Fischer send me this problem				*/
  barrier();

  /* tree configuration */
  if (VERBOSE && me == master)
    printf("Building tree configuration ...\n");
  degree = 0;
  uptree = NULL;
  for (i=0; i<MAXT; i++) downtree[i] = NULL;
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

  /* synchronize because treechannels have to be initialized correctly */
  barrier();

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

  /* initialize asynchronous vchannel array */
  for (i=0; i<ASYNCMAX; i++) {
    asyncadress = GET_AVC_ADR(i);
    SET_AVC_PE(asyncadress,NOPE);
    SET_AVC_ID(asyncadress,NOID);
    SET_AVC_INDEX(asyncadress,NOINDEX);
    SET_AVC_NEXTSEND(asyncadress,0);
    SET_AVC_NEXTRECEIVE(asyncadress,0);
    SET_AVC_CARRIERS(asyncadress,i);
    cararray = GET_AVC_CARRIERS(asyncadress);
    for (j=0; j<ASYNCSENDMAX; j++) {
      carrier = &cararray->out_carriers[j];
      SET_CARRIER_RECEIVED(carrier,RECEIVED);
      SET_CARRIER_SIZE(carrier,NOSIZE);
      SET_CARRIER_DATA(carrier,NULL);
    }
    for (j=0; j<ASYNCRECEIVEMAX; j++) {
      carrier = &cararray->in_carriers[j];
      SET_CARRIER_RECEIVED(carrier,RECEIVED);
      SET_CARRIER_SIZE(carrier,NOSIZE);
      SET_CARRIER_DATA(carrier,NULL);
    }
  }

  /* initialize mailbox */
  mailbox.first_read = 0;
  mailbox.first_write = 0;
  for (i=0; i<MAILMAX; i++) {
    mailadr = GET_MAIL_ADRESS(i);
    SET_MAIL_OCCUPIED(mailadr,NOTOCCUPIED);
    SET_MAIL_SOURCE(mailadr,NOPE);
    SET_MAIL_ID(mailadr,NOID);
    SET_MAIL_SIZE(mailadr,NOSIZE);
    for (j=0; j<RAND_MSG_SIZE; j++) mailadr->data[j] = '0';
    SET_MAIL_READABLE(mailadr,NOTREADABLE);
  }

  /* due to race conditions:											  */
  /* ensure that all procs are initialized before one leaves InitPPIF() */
  barrier();

  return(0);
}

void ExitPPIF (void)
{
  /* to avoid race conditions */
  barrier();

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

  if (procs == 1) return(0);

  error = 0;

  if (me!=master)
    if (RecvSync(uptree,data,size)!=size) error = 1;

  for (i=0; i<degree; i++)
    if (SendSync(downtree[i],data,size)!=size) error = 1;

  return(error);
}

int Concentrate (void *data, int size)
{
  if (me!=master)
  {
    if(SendSync(uptree,data,size)!= size) return(1);
  }
  return(0);
}

int GetConcentrate (int slave, void *data, int size)
{
  if (slave<degree)
  {
    if (RecvSync(downtree[slave],data,size)!=size) return(1);
  }
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
  barrier();
  return(0);
}

/****************************************************************************/
/*                                                                          */
/* Synchronous communication                                                */
/*                                                                          */
/****************************************************************************/

VChannelPtr ConnSync (int p, int id)
{
  int i;
  int ident;            /* channel identification  value                        */
  int prev;                     /* processor, which entered critical section before me  */
  SYNCVCHANNEL *adress;                /* VChannel adress                       */
  VChannelPtr vc;

  /* begin of my side of the deal */

  /* look for free VC Index */
  for (i=0; i<SYNCMAX; i++)
  {
    if (GET_SVCA_PE(i) == NOPE) break;
  }

  /* there is no more free channel */
  if (i >= SYNCMAX)
  {
    printf("ConnSync(): SYNCMAX %d for PE %d exeeded\n", SYNCMAX, me);
    vc = NULL;
    return(vc);
  }

  /* initialize virtual channel */
  adress = GET_SVCA_ADR(sync_channels,i);
  SET_SVC_PE(adress,p);
  SET_SVC_ID(adress,id);
  SET_SVC_INDEX(adress,NOINDEX);
  SET_SVC_DATA(adress,NULL);

  /* make channel identification value */
  ident = MAKE_IDENT(id,i);

  /* Get semaphor value for partner's critical data and set it to me */
  prev = shmem_swap((long *) &sync_flag[me], ident, p);

  /* end of my side of the deal */

  /* begin of partner's side of the deal */

  /* wait until partner has reached equal point and writes his index */
  shmem_wait((long *) &sync_flag[p], NOFLAG);

  /* does p want to connect with same channel id */
  if (GET_IDENT_ID(sync_flag[p]) != id)
  {
    printf("ConnSync(): different channel ids, me=%d myid=%d he=%d hisid=%d\n", me, id, p, GET_IDENT_ID(sync_flag[p]));

    vc = NULL;

    /* exit due to error */
    return(vc);
  }

  /* set channel to corresponding channel index */
  SET_SVC_INDEX(adress,GET_IDENT_INDEX(sync_flag[p]));

  /* reset sync_flag to NOINDEX */
  sync_flag[p] = NOFLAG;

  /* due to race conditions:                   */
  /* tell partner that I am ready and          */
  /* wait until partner has also done his work */
  shmem_swap((long *) &sync_flag[p], me, p);
  shmem_wait((long *) &sync_flag[me], NOFLAG);

  /* reset this flag */
  sync_flag[me] = NOFLAG;

  /* end of partner's side of the deal */

  /* convert to VChannelPtr */
  vc = (VChannelPtr) adress;

  return(vc);
}

int DiscSync (VChannelPtr svc)
{
  int p,index;
  void *adress;
  SYNCVCHANNEL *vc;

  /* convert general pointer */
  vc = (SYNCVCHANNEL *) svc;

  /* store partner's data */
  p = GET_PE(vc);
  index = GET_INDEX(vc);
  adress = GET_SVC_SRECVADR(vc);

  /* recycle virtual channel */
  SET_SVC_PE(vc,NOPE);
  SET_SVC_ID(vc,NOID);
  SET_SVC_INDEX(vc,NOINDEX);
  SET_SVC_DATA(vc,NULL);

  /* reset partner's received flag */
  shmem_swap(adress, RECEIVED, p);

  /* wait until partner sets received_flag */
  shmem_wait((long *) GET_SVC_RECVADR(vc), NOTRECEIVED);

  /* reset received flag */
  SET_SVC_RECEIVED(vc,NOTRECEIVED);

  return(0);
}

int SendSync (VChannelPtr svc, void *data, int size)
{
  void *dataadr;                 /* pointer to data entry */
  void *sizeadr;                 /* pointer to size entry */
  int proc;                      /* proc num of receiver  */
  int sendsize = size;
  SYNCVCHANNEL *vc;

  /* convert general pointer */
  vc = (SYNCVCHANNEL *) svc;

  /* align size to PACKSIZE */
  if (sendsize%PACKSIZE!=0)
    sendsize = sendsize + (PACKSIZE-sendsize%PACKSIZE);

  /* is channel proper initialized */
  if (GET_PE(vc) == NOPE || GET_ID(vc) == NOID || GET_INDEX(vc) == NOINDEX)
  {
    printf("SendSync(): channel not proper initialized me=%d he=%d id=%d index=%d\n", me, GET_PE(vc), GET_ID(vc), GET_INDEX(vc));
    return(-1);
  }

  /* size has to be a multiple of PACKSIZE */
  if (size%PACKSIZE != 0 && PACKSIZE>ALIGNMENT)
  {
    printf("SendSync(): ERROR: size of data not a multiple of PACKSIZE=%d and ALIGNMENT=%d me=%d he=%d size=%d sizerest=%d\n", PACKSIZE, ALIGNMENT, me, GET_PE(vc), size, size%PACKSIZE);
    return(-1);
  }

  /* determine adress to store data pointer at receiver */
  dataadr = GET_SVC_DATAADR(vc);

  /* deterine adress to store size of data  at receiver */
  sizeadr = GET_SVC_SIZEADR(vc);

  /* determine receiver */
  proc = GET_PE(vc);

  /* put data pointer to receiver */
  shmem_swap(dataadr, (long) data, proc);

  /* put size to receiver */
  shmem_swap(sizeadr, sendsize, proc);

  /* output the send operation */
  if (SYNCVERBOSE)
  {
    printf("%4d: SendSync(): dst=%2d id=%5d data=%8x sendsize=%7d\n",me,GET_PE(vc),GET_ID(vc),data,sendsize);
  }

  /* wait until he has gotten the data */
  shmem_wait((long *) GET_SVC_RECVADR(vc), NOTRECEIVED);

  /* reset this flag for next send */
  SET_SVC_RECEIVED(vc,NOTRECEIVED);

  return(size);
}

int RecvSync (VChannelPtr svc, void *data, int size)
{
  int recvsize = size;
  SYNCVCHANNEL *vc;

  /* convert general pointer */
  vc = (SYNCVCHANNEL *) svc;

  /* align size to PACKSIZE */
  if (recvsize%PACKSIZE!=0)
    recvsize = recvsize + (PACKSIZE-recvsize%PACKSIZE);

  /* is channel proper initialized */
  if (GET_PE(vc) == NOPE || GET_ID(vc) == NOID || GET_INDEX(vc) == NOINDEX)
  {
    printf("RecvSync(): channel not proper initialized me=%d he=%d id=%d index=%d\n", me, GET_PE(vc), GET_ID(vc), GET_INDEX(vc));
    return(-1);
  }

  /* size has to be a multiple of PACKSIZE */
  if (size%PACKSIZE != 0 && PACKSIZE>ALIGNMENT)
  {
    printf("RecvSync(): size of data not a multiple of PACKSIZE=%d and ALIGNMENT=%d me=%d he=%d size=%d sizerest=%d\n", PACKSIZE, ALIGNMENT, me, GET_PE(vc), size, size%PACKSIZE);
    return(-1);
  }

  /* wait until sender has written data pointer and size to me */
  shmem_wait((long *) &GET_SIZE(vc), NOSIZE);

  /* compare size of data blocks */
  if (GET_SIZE(vc) != recvsize)
  {
    printf("RecvSync(): different block size, me=%d mysize=%d he=%d hissize=%d\n", me, recvsize, GET_PE(vc), GET_SIZE(vc));

    /* exit due to error */
    return(-1);
  }

        #ifdef TIME
  timebegin = CurrentTime();
        #endif

  /* select packagesize */
        #ifdef PCKSIZE4BYTE
  /* get data from sender */
  shmem_get32(data, GET_SVC_DATA(vc), recvsize/4, GET_PE(vc));
        #else
  /* get data from sender */
  shmem_get(data, GET_SVC_DATA(vc), recvsize/8, GET_PE(vc));
        #endif

        #ifdef TIME
  timeend = CurrentTime();
  printf("PE %d:shmem_get() in RecvSync() sent MAXDATA=%d bytes in time=%f [s]\n",me,recvsize,timeend-timebegin);
        #endif

  /* print the receive parameters */
  if (SYNCVERBOSE)
  {
    printf("%4d: RecvSync(): src=%2d id=%5d data=%8x recvsize=%7d\n",me,GET_PE(vc),GET_ID(vc),data,recvsize);
  }

  /* reset data pointer and size for next receive */
  SET_SVC_DATA(vc,NULL);
  SET_SVC_SIZE(vc,NOSIZE);

  /* now tell sender that I 've reveived the data */
  shmem_swap((long *) GET_SVC_SRECVADR(vc), RECEIVED, GET_PE(vc));

  return(size);
}


/****************************************************************************/
/*                                                                          */
/* Asynchronous communication                                               */
/*                                                                          */
/****************************************************************************/

VChannelPtr ConnASync (int p, int id)
{
  int i,j;
  ASYNCVCCARRIER *cararray;
  CARRIER *carrier;
  ASYNCVCHANNEL *avc_ptr;         /* AVChannel pointer */
  VChannelPtr vc;                 /* VChannel pointer  */

  /* don't connect to me */
  if (p == me || p >= procs)
  {
    printf("PE %d,ConnASync() WARNING do not connect to p=%d id=%d prev=%d\n",me,p,id,*((&p)-1));
    vc = NULL;
    return(vc);
  }


  /* look for free VC Index */
  for (i=0; i<ASYNCMAX; i++)
  {
    if ((GET_AVC_PE(async_channels,i) == NOPE) &&
        (GET_AVC_ID(async_channels,i) == NOID) ) break;
  }

  if (i >= ASYNCMAX)
  {
    printf("ConnASync(): ASYNCMAX=%d for PE=%d exeeded\n", ASYNCMAX, me);
    vc = NULL;
    return(vc);
  }

  /* initialize virtual channel */
  avc_ptr = GET_AVC_ADR(i);
  SET_AVC_INDEX(avc_ptr,NOINDEX);
  SET_AVC_PE(avc_ptr,p);
  SET_AVC_ID(avc_ptr,id);

  vc = (VChannelPtr) avc_ptr;

  /* test whether all carriers are empty */
  cararray = GET_AVC_CARRIERS(avc_ptr);
  for (j=0; j<ASYNCSENDMAX; j++) {
    carrier = &cararray->out_carriers[j];
    if (GET_CARRIER_RECEIVED(carrier) != RECEIVED ||
        GET_CARRIER_SIZE(carrier) != NOSIZE       ||
        GET_CARRIER_DATA(carrier) != NULL)
    {
      printf("PE %d,DiscASync(): sendcarrier %x NOT EMPTY received=%d, data=%x, size=%d\n", me, carrier, GET_CARRIER_RECEIVED(carrier), GET_CARRIER_DATA(carrier), GET_CARRIER_SIZE(carrier));
    }
  }
  for (j=0; j<ASYNCRECEIVEMAX; j++) {
    carrier = &cararray->in_carriers[j];
    if (GET_CARRIER_RECEIVED(carrier) != RECEIVED ||
        GET_CARRIER_SIZE(carrier) != NOSIZE       ||
        GET_CARRIER_DATA(carrier) != NULL)
    {
      printf("PE %d,DiscASync(): receivecarrier %x NOT EMPTY received=%d, data=%x, size=%d\n", me, carrier, GET_CARRIER_RECEIVED(carrier), GET_CARRIER_DATA(carrier), GET_CARRIER_SIZE(carrier));
    }
  }

  return(vc);
}

int InfoAConn (VChannelPtr avc)
{
  int i,j,index;
  int *index_ptr;
  ASYNCVCHANNEL *vc;

  /* convert general pointer */
  vc = (ASYNCVCHANNEL *) avc;

  /* test if partner has completed connection */
  if (GET_INDEX(vc) != NOINDEX) return(1);

  /* fetch async_channel structure from partner */
  shmem_get((long *) remote_channels, (long *) async_channels, ASYNCVCHANNELSIZE/8, GET_PE(vc));

  /* search for corresponding entry */
  for (j=0; j<ASYNCMAX; j++)
  {
    /* test whether this is the searched channel */
    if ( (GET_AVC_PE(remote_channels,j) == me) &&
         (GET_AVC_ID(remote_channels,j) == GET_ID(vc))   )
      break;
  }

  /* has partner already done his ConnASync() */
  if (j < ASYNCMAX)
  {
    /* partner has already established channel with current id */
    /* set my index to value of corresponding entry */
    SET_AVC_INDEX(vc,j);

    /* calculate my index in channelarray */
    i = vc-async_channels;

    /* set his index to my corresponding entry */
    shmem_swap((long *) GET_AVC_INDEXADR(async_channels,j), i, GET_PE(vc));

    /* wait until data has arrived */
    _remote_write_barrier();

    /* connection complete */
    return(1);
  }
  else
  {
    /* connection still in progress */
    return(0);
  }
  return(-1);
}

int DiscASync (VChannelPtr avc)
{
  int index,j;
  int pe;
  ASYNCVCCARRIER *cararray;
  CARRIER *carrier;
  ASYNCVCHANNEL *vc;

  /* convert general pointer */
  vc = (ASYNCVCHANNEL *) avc;

  /* check whether channel is proper initialized */
  /* id may be already reset by partner          */
  if (GET_PE(vc) == NOPE || GET_INDEX(vc) == NOINDEX)
  {
    printf("PE %d:DiscASync(): channel not proper initialized me=%d he=%d id=%d index=%d\n", me, me, GET_PE(vc), GET_ID(vc), GET_INDEX(vc));
    return(-1);
  }

  /* test whether all carriers are empty */
  cararray = GET_AVC_CARRIERS(vc);
  for (j=0; j<ASYNCSENDMAX; j++) {
    carrier = &cararray->out_carriers[j];
    if (GET_CARRIER_RECEIVED(carrier) != RECEIVED ||
        GET_CARRIER_SIZE(carrier) != NOSIZE       ||
        GET_CARRIER_DATA(carrier) != NULL)
    {
      printf("PE %d,DiscASync(): sendcarrier %x NOT EMPTY received=%d, data=%x, size=%d\n", me, carrier, GET_CARRIER_RECEIVED(carrier), GET_CARRIER_DATA(carrier), GET_CARRIER_SIZE(carrier));
    }
  }
  for (j=0; j<ASYNCRECEIVEMAX; j++) {
    carrier = &cararray->in_carriers[j];
    if (GET_CARRIER_RECEIVED(carrier) != RECEIVED ||
        GET_CARRIER_SIZE(carrier) != NOSIZE       ||
        GET_CARRIER_DATA(carrier) != NULL)
    {
      printf("PE %d,DiscASync(): receivecarrier %x NOT EMPTY received=%d, data=%x, size=%d\n", me, carrier, GET_CARRIER_RECEIVED(carrier), GET_CARRIER_DATA(carrier), GET_CARRIER_SIZE(carrier));
    }
  }

  /* recycle VChannel */
  index = GET_AVC_INDEX(vc);
  pe = GET_PE(vc);
  SET_AVC_INDEX(vc,NOINDEX);
  SET_AVC_PE(vc,NOPE);
  SET_AVC_NEXTSEND(vc,0);
  SET_AVC_NEXTRECEIVE(vc,0);
  shmem_swap((long *) GET_AVC_IDADR(async_channels,index),NOID,pe);

  /* wait until data has arrived */
  _remote_write_barrier();

  return(0);
}

int InfoADisc (VChannelPtr avc)
{
  ASYNCVCHANNEL *vc;

  /* convert general pointer */
  vc = (ASYNCVCHANNEL *) avc;

  /* if partner has also done his DiscASync() */
  /* then VChannel id is reset                */
  if (GET_ID(vc) == NOID)
  {
    /* disconnection complete */
    return(1);
  }
  else
  {
    /* disconnection still in progress */
    return(0);
  }
  return(-1);
}


msgid SendASync (VChannelPtr avc, void *data, int size, int *error)
{
  int index;
  CARRIER *carrier;
  ASYNCVCHANNEL *vc;

  /* convert general pointer */
  vc = (ASYNCVCHANNEL *) avc;

  /* is channel proper initialized */
  if (GET_PE(vc) == NOPE || GET_ID(vc) == NOID || GET_INDEX(vc) == NOINDEX)
  {
    printf("SendASync(): channel not proper initialized me=%d he=%d id=%d index=%d\n", me, GET_PE(vc), GET_ID(vc), GET_INDEX(vc));
    *error = 1;
    return(NULL);
  }

  /* size has to be a multiple of PACKSIZE */
  if (size%PACKSIZE != 0)
  {
    printf("SendASync(): size of data not a multiple of PACKSIZE=%d me=%d he=%d size=%d sizerest=%d\n", PACKSIZE, me, GET_PE(vc), size, size%PACKSIZE);
    *error = 1;;
    return(NULL);
  }

  /* output the send operation */
  if (ASYNCVERBOSE)
  {
    printf("%4d: SendASync(): dst=%2d id=%5d data=%8x size=%7d\n",me,GET_PE(vc),GET_ID(vc),data,size);
  }

  /* get carrier index */
  index = GET_AVC_NEXTSEND(vc);

  /* get pointer to carrier */
  carrier = GET_SCARRIER_ADR(vc,index);

  /* is carrier really empty */
  if (GET_CARRIER_RECEIVED(carrier) != RECEIVED ||
      GET_CARRIER_DATA(carrier) != NULL ||
      GET_CARRIER_SIZE(carrier) != NOSIZE )
  {
    printf("SendASync(): next send entry=%d of AVChannel  not empty on PE %d\n", index, me);
    *error = 1;
    return(NULL);
  }

  /* initialize carrier */
  SET_CARRIER_RECEIVED(carrier,NOTRECEIVED);
  SET_CARRIER_SIZE(carrier,size);
  SET_CARRIER_DATA(carrier,data);

  /* increment next entry index */
  SET_AVC_NEXTSEND(vc,(GET_AVC_NEXTSEND(vc)+1)%ASYNCSENDMAX);

  /* no error occured */
  *error = 0;

  return(carrier);
}


int InfoASend (VChannelPtr vc, msgid m)
{

  /* check if receiver got the data */
  if (GET_CARRIER_RECEIVED(m) == RECEIVED)
  {
    /* transfer complete */
    return(1);
  }
  else
  {
    /* transfer in progress */
    return(0);
  }
  return(-1);
}


msgid RecvASync (VChannelPtr avc, void *data, int size, int *error)
{
  int index;
  CARRIER *carrier;
  ASYNCVCHANNEL *vc;

  /* convert general pointer */
  vc = (ASYNCVCHANNEL *) avc;

  /* is channel proper initialized */
  if (GET_PE(vc) == NOPE || GET_ID(vc) == NOID || GET_INDEX(vc) == NOINDEX)
  {
    printf("RecvASync(): channel not proper initialized me=%d he=%d id=%d index=%d\n", me, GET_PE(vc), GET_ID(vc), GET_INDEX(vc));
    *error = 1;
    return(NULL);
  }

  /* print the receive parameters */
  if (ASYNCVERBOSE)
  {
    printf("%4d: RecvASync(): src=%2d id=%5d data=%8x size=%7d\n",me,GET_PE(vc),GET_ID(vc),data,size);
  }

  /* size has to be a multiple of PACKSIZE */
  if (size%PACKSIZE != 0)
  {
    printf("RecvASync(): size of data not a multiple of PACKSIZE=%d me=%d he=%d size=%d sizerest=%d\n", PACKSIZE, me, GET_PE(vc), size, size%PACKSIZE);
    *error = 1;
    return(NULL);
  }

  /* determine next receive entry */
  index = GET_AVC_NEXTRECEIVE(vc);

  /* get pointer to carrier */
  carrier = GET_RCARRIER_ADR(vc,index);

  /* is carrier really empty */
  if (GET_CARRIER_RECEIVED(carrier) != RECEIVED)
  {
    printf("RecvASync(): next receive entry=%d of AVChannel not empty on PE %d\n", index, me);
    *error = 1;
    return(NULL);
  }

  /* initialize carrier */
  SET_CARRIER_SIZE(carrier,size);
  SET_CARRIER_DATA(carrier,data);

  /* increment messages to receive */
  SET_AVC_NEXTRECEIVE(vc,(GET_AVC_NEXTRECEIVE(vc)+1)%ASYNCRECEIVEMAX);

  /* no error occured */
  *error = 0;

  return(carrier);
}


int InfoARecv (VChannelPtr avc, msgid m)
{
  int index;
  ASYNCVCHANNEL *sendvc;
  CARRIER *scar_ptr;
  ASYNCVCHANNEL *vc;

  /* convert general pointer */
  vc = (ASYNCVCHANNEL *) avc;

  /* compute corresponding VChannel adress */
  sendvc = GET_AVC_ADR(GET_INDEX(vc));

  /* compute corresponding send carrier index */
  index = m-GET_RCARRIER_ADR(vc,0);

  /* compute address of correspondig send carrier adress */
  scar_ptr = GET_SCARRIER_ADR(sendvc,index);

  /* get corresponding send carrier */
  shmem_get((long *) carrier, (long *) scar_ptr, sizeof(CARRIER)/8, GET_PE(vc));

  /* check whether partner has done his ASyncSend() */
  if (GET_CARRIER_DATA(carrier) != NULL && GET_CARRIER_SIZE(carrier) != NOSIZE)
  {
    /* compare size of datablock, etc. */
    if (GET_CARRIER_SIZE(carrier) != GET_CARRIER_SIZE(m) ||
        GET_PE(vc)>procs-1                               ||
        GET_PE(vc)<0                                     ||
        GET_CARRIER_SIZE(m) <= 0                         ||
        GET_CARRIER_DATA(m) == NULL                      ||
        GET_CARRIER_DATA(carrier) == NULL                    ||
        GET_CARRIER_SIZE(m)%PACKSIZE != 0)
    {
      printf("PE %d:InfoARecv(): he=%d VC=%x vcindex=%d msgid=%x index=%d ERROR diff size: mysize=%d hissize=%d mydata=%x hisdata=%x\n", me,  GET_PE(vc), vc, vc-async_channels, m, index, GET_CARRIER_SIZE(m), GET_CARRIER_SIZE(carrier),GET_CARRIER_DATA(m), GET_CARRIER_DATA(carrier));
      return(-1);
    }

                #ifdef PCKSIZE4BYTE
    /* get data from sender */
    shmem_get32((short *) GET_CARRIER_DATA(m), (short *) GET_CARRIER_DATA(carrier), (int) GET_CARRIER_SIZE(m)/4, GET_PE(vc));
                #else
    /* get data from sender */
    shmem_get((long *) GET_CARRIER_DATA(m), (long *) GET_CARRIER_DATA(carrier), GET_CARRIER_SIZE(m)/8, GET_PE(vc));
                #endif

    /* recylce carrier at my side */
    SET_CARRIER_DATA(m,NULL);
    SET_CARRIER_SIZE(m,NOSIZE);
    SET_CARRIER_RECEIVED(m,RECEIVED);

    /* recycle allocated carrier */
    SET_CARRIER_DATA(carrier,NULL);
    SET_CARRIER_SIZE(carrier,NOSIZE);
    SET_CARRIER_RECEIVED(carrier,NOTRECEIVED);

    /* recycle carrier at sender's side */
    shmem_put((long *) scar_ptr, (long *) carrier, sizeof(CARRIER)/8, GET_PE(vc));

    /* set received flag at sender */
    shmem_swap((long *) &GET_CARRIER_RECEIVED(scar_ptr), RECEIVED, GET_PE(vc));

    /* wait until data has arrived */
    _remote_write_barrier();

    /* reset this flag */
    SET_CARRIER_RECEIVED(carrier,RECEIVED);

    /* data transfer complete */
    return(1);

  }
  else
  {
    /* recycle allocated carrier */
    /*SET_CARRIER_DATA(carrier,NULL);
       SET_CARRIER_SIZE(carrier,NOSIZE);
       SET_CARRIER_RECEIVED(carrier,RECEIVED); */

    /* compare size of datablock, etc. */
    if (ASYNCVERBOSE &&
        (GET_PE(vc)>procs-1                   ||
         GET_PE(vc)<0                         ||
         GET_CARRIER_SIZE(m) <= 0             ||
         GET_CARRIER_DATA(m) == NULL          ||
         GET_CARRIER_DATA(carrier) != NULL        ||
         GET_CARRIER_SIZE(carrier) != NOSIZE  ||
         GET_CARRIER_SIZE(m)%PACKSIZE != 0))
    {
      printf("PE %d:InfoARecv(): he=%d VC=%x vcindex=%d msgid=%x index=%d WARNING carrier inconsistent: mysize=%d hissize=%d mydata=%x hisdata=%x\n", me,  GET_PE(vc), vc, vc-async_channels, m, index, GET_CARRIER_SIZE(m), GET_CARRIER_SIZE(carrier), GET_CARRIER_DATA(m), GET_CARRIER_DATA(carrier));
    }


    /* data transfer still in progress */
    return(0);
  }
  return(-1);
}


/****************************************************************************/
/*                                                                          */
/* Random communication                                                     */
/*                                                                          */
/****************************************************************************/


int SendMail (int destId, int reqId, void *data, int size)
{
  int first_write;
  int i,index;
  char *data_byte;
  char *maildata_byte;

  /* get first write entry of destination mailbox */
  shmem_get((long *) &first_write, (long *) &mailbox.first_write, 1, destId);

  /* search for entry which is NOW free (result of concurrent access) */
  index = 0;
  while ((shmem_swap((long *) GET_MAIL_ADRESS(first_write), OCCUPIED, destId)
          == OCCUPIED) && (index<MAILMAX)) {
    first_write = NEXT_MAIL_INDEX(first_write);
    index++;
  }

  /* I 've lost during concurrent acess */
  if (index >= MAILMAX)
  {
    printf("During concurrent access: Maximum of mails %d for PE %d exeeded\n", MAILMAX, destId);
    return(1);
  }

  /* now free is my entry to use for mail storage                */
  /* due to race conditions mail ordering might NOT be preserved */

  /* fill out mail to receiver */
  SET_MAIL_OCCUPIED(newmail,OCCUPIED);
  SET_MAIL_SOURCE(newmail,me);
  SET_MAIL_ID(newmail,reqId);
  SET_MAIL_SIZE(newmail,size);

  /* copy data to data of mail */
  data_byte = (char *) data;
  maildata_byte = (char *) &GET_MAIL_DATA(newmail);
  for (i=0; i<size; i++)
    *(maildata_byte++) = *(data_byte++);

  SET_MAIL_READABLE(newmail,NOTREADABLE);

        #ifdef PCKSIZE4BYTE
  /* send message */
  shmem_put32((short *) GET_MAIL_ADRESS(first_write), (short *) newmail, MAILSIZE/4, destId);
        #else
  /* send message */
  shmem_put((long *) GET_MAIL_ADRESS(first_write), (long *) newmail, MAILSIZE/8, destId);
        #endif

  /* set readable flag */
  shmem_swap((long *) &GET_MAIL_READABLE(first_write), READABLE, destId);

  /* recycle mailbox entry */
  SET_MAIL_SOURCE(newmail,NOPE);
  SET_MAIL_ID(newmail,NOID);
  SET_MAIL_SIZE(newmail,NOSIZE);
  for (i=0; i<GET_MAIL_SIZE(newmail); i++)
    *maildata_byte++ = NODATA;
  SET_MAIL_READABLE(newmail,NOTREADABLE);
  SET_MAIL_OCCUPIED(newmail,NOTOCCUPIED);

  return(0);
}

int GetMail (int *sourceId, int *reqId, void *data, int *size)
{
  int first_read;
  int cnt;
  int i,j;
  int first_readable;
  int first_writeable;
  char *data_byte;
  char *maildata_byte;
  MAIL *first_ptr;

  /* initialize to undefined */
  first_readable = NOFIRST;
  first_writeable = NOFIRST;

  /* no messages in mailbox */
  cnt = 0;

  /* determine index to read next */
  first_read = mailbox.first_read;

  /* look for first message and count messages */
  for (i=first_read,j=0; j<MAILMAX; i=NEXT_MAIL_INDEX(i),j++) {
    if (GET_MAIL_READABLE(i) == READABLE)
    {
      if (first_readable==NOFIRST)
      {
        first_readable = i;
      }
      cnt++;
    }
    if (GET_MAIL_OCCUPIED(GET_MAIL_ADRESS(i)) == NOTOCCUPIED &&
        first_readable != NOFIRST &&
        first_writeable == NOFIRST)
    {
      first_writeable = i;
    }
  }

  /* mailbox empty leave everything unchanged */
  if (cnt == 0)
  {
    return(0);
  }

  /* mail ordering because of race conditions not preserved */
  if (first_readable != first_read)
  {
    printf("GetMail(): Mail ordering for PE %d is not preserved: first_read=%d, first_readable=%d\n", me, first_read, first_readable);
  }

  /* extract mail entry */
  first_ptr = GET_MAIL_ADRESS(first_readable);
  *sourceId = GET_MAIL_SOURCE(first_ptr);
  *reqId = GET_MAIL_ID(first_ptr);
  *size = GET_MAIL_SIZE(first_ptr);

  /* copy data of mail to data */
  data_byte = (char *) data;
  maildata_byte = (char *) GET_MAIL_DATA(first_ptr);
  for (i=0; i<GET_MAIL_SIZE(first_ptr); i++)
    *(data_byte++) = *(maildata_byte++);

  /* recycle mailbox entry */
  SET_MAIL_SOURCE(first_ptr,NOPE);
  SET_MAIL_ID(first_ptr,NOID);
  SET_MAIL_SIZE(first_ptr,NOSIZE);
  for (i=0; i<GET_MAIL_SIZE(first_ptr); i++)
    *maildata_byte++ = NODATA;
  SET_MAIL_READABLE(first_ptr,NOTREADABLE);
  SET_MAIL_OCCUPIED(first_ptr,NOTOCCUPIED);

  /* set mailbox status */
  mailbox.first_read = NEXT_MAIL_INDEX(first_readable);
  if (first_writeable!=NOFIRST)
    mailbox.first_write = first_writeable;
  else
    mailbox.first_write = mailbox.first_read;

  return(cnt);
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
  return(((double)(rtclock())/((double)clocks)));
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
