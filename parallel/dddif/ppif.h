// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ppif.h                                                        */
/*                                                                          */
/* Purpose:   header file for parallel processor interface                  */
/*                                                                          */
/* Author:    Peter Bastian / Klaus Birken                                  */
/*                                                                          */
/* History:   17 Aug 1992 begin                                             */
/*            14 Sep 1993 pass argc, argv from main to InitPPIF             */
/*            16 Sep 1993 async send/recv return msgid now                  */
/*            951106 kb  changed parameters for InitPPIF()                  */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __PPIF__
#define __PPIF__
#ifdef __cplusplus
extern "C" {
#endif


/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef int *VChannelPtr;   /* dummy definition, any pointer type is ok     */
typedef int msgid;          /* type of return value for async comm !32bit!  */

enum directions {north,east,south,west,up,down};

#define RAND_MSG_SIZE 128   /* max size of random messages                                      */

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/* id's */
extern int me;                     /* my processor id                       */
extern int master;                 /* id of master processor                */
extern int procs;                  /* number of processors in the network   */

/* 3D array structure */
extern int arrayid;                                /* position in compact form 8 bits each  */
extern int MyX,MyY,MyZ;            /* 3D array coordinates                  */
extern int DimX,DimY,DimZ;         /* 3D array dimensions, may be 1 !       */
extern VChannelPtr nn[6];          /* nearest neighbors in 3D array         */

/* Tree structure */
extern int degree;                 /* degree of downtree nodes              */
extern VChannelPtr uptree;         /* channel uptree                        */
extern VChannelPtr downtree[];     /* channels downtree (may be empty)      */
extern int slvcnt[];                   /* number of processors in subtree       */


/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/

#define APOS_TO_AID(x,y,z)  ((z<<16)|(y<<8)|x)   /* array pos to compact for*/
#define XPOS(aid)           (aid&0xFF)           /* xpos from compact form  */
#define YPOS(aid)           ((aid&0xFF00)>>8)    /* ypos from compact form  */
#define ZPOS(aid)           ((aid&0xFF0000)>>16) /* zpos from compact form  */

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/* initialization & shutdown */
int         InitPPIF         (int *argcp, char ***argvp);
int         ExitPPIF         (void);

/* tree oriented functions */
int         Broadcast        (void *data, int size);
int         Concentrate      (void *data, int size);
int         GetConcentrate   (int slave, void *data, int size);
int         Spread           (int slave, void *data, int size);
int         GetSpread        (void *data, int size);
int         Synchronize      (void);

/* synchronous communication */
VChannelPtr ConnSync         (int p, int id);
int         DiscSync         (VChannelPtr vc);
int         SendSync         (VChannelPtr vc, void *data, int size);
int         RecvSync         (VChannelPtr vc, void *data, int size);

/* asynchronous communication */
VChannelPtr ConnASync        (int p, int id);
int         DiscASync        (VChannelPtr vc);
msgid       SendASync        (VChannelPtr vc, void *data, int size, int *error);
msgid       RecvASync        (VChannelPtr vc, void *data, int size, int *error);
int         InfoAConn        (VChannelPtr vc);
int         InfoADisc        (VChannelPtr vc);
int         InfoASend        (VChannelPtr vc, msgid m);
int         InfoARecv        (VChannelPtr vc, msgid m);

/* random communication */
int             SendMail                 (int destId, int reqId, void *data, int size);
int             GetMail                  (int *sourceId, int *reqId, void *data, int *size);

/* miscellaneous */
int         UsedSpace        (void);
void        PrintHostMessage (char *s);
double      CurrentTime      (void);
int         Distance             (int p, int q);
int         aid_to_pid       (int x, int y, int z);
int         pid_to_aid       (int p);


#ifdef __cplusplus
}
#endif
#endif
