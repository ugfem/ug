// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      lowcomm.h                                                     */
/*                                                                          */
/* Purpose:   lowlevel communication layer                                  */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   960715 kb  begin                                              */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/* RCS_ID
   $Header: /hosts/dom/cvs/ddd/src/dddi.h,v 1.8 1997/07/24 15:49:05 birken Exp
   $
 */


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __DDD_LOWCOMM_H__
#define __DDD_LOWCOMM_H__



/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define EXCEPTION_LOWCOMM_CONNECT  -10

/* lowcomm users should use exceptions EXCEPTION_LOWCOMM_USER or lower */
#define EXCEPTION_LOWCOMM_USER     -100


/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/


typedef unsigned long ULONG;


typedef void *LC_MSGHANDLE;  /* handle for actual message (send OR recv side*/
typedef void *LC_MSGTYPE;    /* type of message (on send AND recv side) */
typedef int LC_MSGCOMP;      /* component of message (dto) */


/* function pointer types for alloc and free */
typedef void * (*AllocFunc)(size_t);
typedef void (*FreeFunc)(void *);


/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/


/* lowcomm.c */
void  LC_Init (AllocFunc,FreeFunc);
void  LC_Exit (void);

void  LC_SetMemMgrSend (AllocFunc,FreeFunc);
void  LC_SetMemMgrRecv (AllocFunc,FreeFunc);
void  LC_SetMemMgrDefault (void);


LC_MSGTYPE LC_NewMsgType (char *);
LC_MSGCOMP LC_NewMsgTable (char *, LC_MSGTYPE, size_t);
LC_MSGCOMP LC_NewMsgChunk (char *, LC_MSGTYPE);

void       LC_MsgSend (LC_MSGHANDLE);

int           LC_Connect (LC_MSGTYPE);
int           LC_Abort (int);
LC_MSGHANDLE *LC_Communicate (void);
void          LC_Cleanup (void);



LC_MSGHANDLE LC_NewSendMsg (LC_MSGTYPE, DDD_PROC);
ULONG    LC_GetTableLen (LC_MSGHANDLE, LC_MSGCOMP);
void *   LC_GetPtr (LC_MSGHANDLE, LC_MSGCOMP);
DDD_PROC LC_MsgGetProc (LC_MSGHANDLE);

size_t   LC_MsgPrepareSend (LC_MSGHANDLE);
size_t   LC_MsgFreeze (LC_MSGHANDLE);
int      LC_MsgAlloc(LC_MSGHANDLE);

void     LC_SetTableLen (LC_MSGHANDLE, LC_MSGCOMP, ULONG);
void     LC_SetTableSize (LC_MSGHANDLE, LC_MSGCOMP, ULONG);
void     LC_SetChunkSize (LC_MSGHANDLE, LC_MSGCOMP, size_t);

size_t   LC_GetBufferSize (LC_MSGHANDLE);


void LC_PrintSendMsgs (void);
void LC_PrintRecvMsgs (void);


#endif
